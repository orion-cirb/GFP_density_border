package GFP_density_border_Tools;

import GFP_density_border_Tools.Cellpose.CellposeTaskSettings;
import GFP_density_border_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;

import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import sc.fiji.skeletonize3D.Skeletonize3D_;
import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Voxel3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureFeret;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/**
 * Tools for bacteria / DNA detection (Omnipose) and border-fragment GFP/DAPI intensity analysis.
 *
 * @author Orion-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));

    public Calibration cal = new Calibration();
    private double pixelSurf = 0;
    String[] channelsName = {"Bacteria: ", "DNA: ", "GFP: "};

    // Omnipose
    private String omniposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home") + "\\miniconda3\\envs\\omnipose\\" :
            "/opt/miniconda3/envs/omnipose/";
    private String omniposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home") + "\\.cellpose\\models\\" :
            System.getProperty("user.home") + "/.cellpose/models/";
    public String omniposeBactModel = "bact_phase_omnitorch_0";
    public String omniposeDnaModel = "bact_fluor_omnitorch_0";

    private int omniposeDiameter = 18;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0;
    private boolean useGpu = true;

    // Bacteria size filtering
    public double minBactSurface = 0.4;
    public double maxBactSurface = 20;

    // DNA size filtering
    public double minDnaSurface = 0.4;
    public double maxDnaSurface = 20;

    // Border fragments parameters
    public double borderWidth = 1;
    public int nFragments = 10;

    // Total DNA (or shared-DNA) volume associated to each bacterium label, used for volume ratio computation
    private Map<Integer, Double> dnaTotalBactVol = new HashMap<>();


    /**
     * Print a message to both the ImageJ console and the status bar.
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }


    /**
     * Check that the 3D ImageJ Suite is installed (required dependency).
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }


    /**
     * Flush and close an image to free memory.
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }


    /**
     * Detect the image file extension used in the input folder (nd, czi, lif, ics, ics2, lsm, tif, tiff).
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd":
                case "czi":
                case "lif":
                case "ics":
                case "ics2":
                case "lsm":
                case "tif":
                case "tiff":
                    ext = fileExt;
                    break;
            }
        }
        return (ext);
    }


    /**
     * List all image files in a folder matching the given extension.
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList<>();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return (images);
    }


    /**
     * Read XY/Z pixel calibration from OME metadata.
     */
    public void findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }


    /**
     * Retrieve channel names, using a format-specific strategy since metadata fields vary by format.
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt = FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd":
            case "nd2":
                for (int n = 0; n < chs; n++)
                    channels[n] = (meta.getChannelID(0, n) == null) ? Integer.toString(n) : meta.getChannelName(0, n);
                break;
            case "lif":
                for (int n = 0; n < chs; n++)
                    channels[n] = (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                            ? Integer.toString(n) : meta.getChannelName(0, n);
                break;
            case "czi":
                for (int n = 0; n < chs; n++)
                    channels[n] = (meta.getChannelID(0, n) == null) ? Integer.toString(n) : meta.getChannelFluor(0, n);
                break;
            case "ics":
            case "ics2":
                for (int n = 0; n < chs; n++)
                    channels[n] = (meta.getChannelID(0, n) == null) ? Integer.toString(n)
                            : meta.getChannelExcitationWavelength(0, n).value().toString();
                break;
            default:
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return (channels);
    }


    /**
     * Show the parameters dialog box (channels, Omnipose settings, size filters,
     * border-fragment parameters, calibration) and read back the user's choices.
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets(0, 160, 0);
        gd.addImage(icon);

        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channelsName) {
            gd.addChoice(ch, channels, channels[index]);
            index++;
        }

        gd.addMessage("Bacteria and DNA detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath);
        gd.addMessage("Object size threshold ", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min bacterium surface (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium surface (µm2): ", maxBactSurface);
        gd.addNumericField("Min DNA surface (µm2): ", minDnaSurface);
        gd.addNumericField("Max DNA surface (µm2): ", maxDnaSurface);

        gd.addMessage("Border fragments", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Border width (µm): ", borderWidth);
        gd.addNumericField("Number of fragments: ", nFragments);

        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.showDialog();

        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if (gd.wasCanceled())
            ch = null;

        omniposeEnvDirPath = gd.getNextString();
        omniposeModelsPath = gd.getNextString();
        minBactSurface = (float) gd.getNextNumber();
        maxBactSurface = (float) gd.getNextNumber();
        minDnaSurface = (float) gd.getNextNumber();
        maxDnaSurface = (float) gd.getNextNumber();
        borderWidth = gd.getNextNumber();
        nFragments = (int) gd.getNextNumber();

        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixelSurf = cal.pixelWidth * cal.pixelWidth;

        return (ch);
    }


    /**
     * Perform a Z-projection of a stack using the given projection method (e.g. average intensity).
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
        return (zproject.getProjection());
    }


    /**
     * Run Omnipose segmentation on an image and return the resulting population of 3D objects,
     * filtered by size (and optionally excluding objects touching image borders).
     */
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact, String model, double min, double max, boolean excludeBorders) {
        ImagePlus imgIn = new Duplicator().run(imgBact);
        imgIn.setCalibration(cal);

        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath + model, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setOmni(true);
        settings.useMxNet(false);
        settings.setCluster(true);
        settings.setCellProbTh(omniposeMaskThreshold);
        settings.setFlowTh(omniposeFlowThreshold);
        settings.useGpu(useGpu);

        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);

        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        if (excludeBorders)
            pop = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(imgOut), false);
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(min / pixelSurf, max / pixelSurf);
        pop.resetLabels();

        flush_close(imgIn);
        return (pop);
    }


    /**
     * Reset per-image accumulator state. Must be called once per new image/timepoint,
     * before processing bacteria and DNA populations.
     */
    public void initImageVariables() {
        dnaTotalBactVol.clear();
    }


    /**
     * Link each DNA object to its parent bacterium (or bacteria) by testing whether the
     * bacterium centroid falls inside the DNA object. Stores the total bacterial volume
     * associated with each DNA label (single bacterium, or sum if shared by several).
     */
    public void dnaBactLink(Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop) {
        if (bactPop.getNbObjects() == 0 || dnaPop.getNbObjects() == 0) return;

        for (Object3DInt dna : dnaPop.getObjects3DInt()) {
            Map<Integer, Double> matchingBacts = new LinkedHashMap<>();

            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                MeasureCentroid bactCenter = new MeasureCentroid(bact);
                if (dna.contains(bactCenter.getCentroidRoundedAsVoxelInt())) {
                    matchingBacts.put((int) bact.getLabel(), new MeasureVolume(bact).getVolumeUnit());
                }
            }

            if (matchingBacts.size() == 1) {
                dna.setIdObject(matchingBacts.keySet().iterator().next());
                dnaTotalBactVol.put((int) dna.getIdObject(), matchingBacts.values().iterator().next());
            } else if (matchingBacts.size() > 1) {
                dna.setIdObject(matchingBacts.keySet().iterator().next());
                double totalVol = matchingBacts.values().stream().mapToDouble(Double::doubleValue).sum();
                dnaTotalBactVol.put((int) dna.getIdObject(), totalVol);
            }
        }
    }

    /**
     * Assign each DNA object to the bacterium whose mask contains the DNA centroid,
     * then discard DNA objects that don't fall inside any bacterium.
     */
    public void dnaBactNumber(Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop) {
        if (bactPop.getNbObjects() != 0 && dnaPop.getNbObjects() != 0) {
            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                for (Object3DInt dna : dnaPop.getObjects3DInt()) {
                    MeasureCentroid dnaCenter = new MeasureCentroid(dna);
                    if (bact.contains(dnaCenter.getCentroidRoundedAsVoxelInt())) {
                        dna.setIdObject(bact.getLabel());
                    }
                }
            }
        }
        dnaPop.getObjects3DInt().removeIf(p -> p.getIdObject() == 0);
    }

    /**
     * Return the sub-population of DNA objects linked to a given bacterium label.
     */
    private Objects3DIntPopulation findDnaInBact(float bactLabel, Objects3DIntPopulation dnaPop) {
        Objects3DIntPopulation dnaBactPop = new Objects3DIntPopulation();
        for (Object3DInt dna : dnaPop.getObjects3DInt()) {
            if (dna.getIdObject() == bactLabel)
                dnaBactPop.addObject(dna);
        }
        dnaBactPop.resetLabels();
        return (dnaBactPop);
    }


    /**
     * Write the header row of the border-fragment ratio results file.
     * One column per fragment, containing the DAPI/GFP total-intensity ratio.
     * Must be called once, before the main image loop.
     */
    public void writeFragmentHeader(BufferedWriter file, int nFragments) throws IOException {
        StringBuilder header = new StringBuilder("Image\tTime\tBacteria");
        for (int f = 1; f <= nFragments; f++) {
            header.append("\tFragment").append(f).append("_border_DAPI_sum");
            header.append("\tFragment").append(f).append("_border_GFP_sum");
            header.append("\tFragment").append(f).append("_border_DAPI_GFP_ratio");
            header.append("\tFragment").append(f).append("_interior_DAPI_sum");
            header.append("\tFragment").append(f).append("_interior_GFP_sum");
            header.append("\tFragment").append(f).append("_interior_DAPI_GFP_ratio");
        }
        header.append("\n");
        file.write(header.toString());
        file.flush();
    }

    /**
     * For each bacterium: compute its extended skeleton (drawn into the returned image),
     * divide its border into equal-size angular fragments around the centroid, compute
     * the DAPI/GFP total-intensity ratio in each fragment, write one result row per
     * bacterium to fragFile, and save a labeled fragment image to disk.
     *
     * @param img1          reference image (phase/bacteria channel), used for dimensions/calibration
     * @param bactPop       segmented bacteria population
     * @param gfpImg        GFP channel image (denominator of the ratio)
     * @param dapiImg       DAPI/DNA channel image (numerator of the ratio)
     * @param borderWidthUm border thickness in µm
     * @param nFragments    number of equal-size border fragments per bacterium
     * @param imgName       image/series name, written in the results file
     * @param time          timepoint index, written in the results file
     * @param fragFile      writer for the fragment results file (header already written)
     * @param outDir        output directory for the fragment label image
     * @return an ImagePlus containing the extended skeletons of all bacteria
     */
    public ImagePlus bactBorderIntensities(ImagePlus img1, Objects3DIntPopulation bactPop,
                                           ImagePlus gfpImg, ImagePlus dapiImg,
                                           double borderWidthUm, int nFragments,
                                           String imgName, int time, BufferedWriter fragFile, String outDir) throws IOException {

        ImagePlus labelImg = ImageHandler.wrap(img1).createSameDimensions().getImagePlus();
        bactPop.drawInImage(ImageHandler.wrap(labelImg));

        ImagePlus skelResult = ImageHandler.wrap(img1).createSameDimensions().getImagePlus();
        ImagePlus fragResult = ImageHandler.wrap(img1).createSameDimensions().getImagePlus();

        int nBact = bactPop.getNbObjects();
        double pixelSizeUm = img1.getCalibration().pixelWidth;
        ImageProcessor gfpIP = gfpImg.getProcessor();
        ImageProcessor dapiIP = dapiImg.getProcessor();
        ImageProcessor labelIP = labelImg.getProcessor();
        double gfpBackground = estimateBackground(gfpIP, labelIP);
        double dapiBackground = estimateBackground(dapiIP, labelIP);
        System.out.println("Estimated background - GFP: " + gfpBackground + ", DAPI: " + dapiBackground);

        for (int i = 1; i <= nBact; i++) {

            ImagePlus bact = new Duplicator().run(labelImg);
            IJ.setThreshold(bact, i, i);
            IJ.run(bact, "Convert to Mask", "background=Dark black");

            ImageProcessor maskIP = bact.getProcessor().duplicate();

            // Skeletonize and extend the skeleton to reach the true bacterium tips
            Skeletonize3D_ skel = new Skeletonize3D_();
            skel.setup("", bact);
            skel.run(bact.getProcessor());
            ImageProcessor skelIP = bact.getProcessor();
            extendSkeletonEndpoints(skelIP, maskIP);

            int w = skelIP.getWidth(), h = skelIP.getHeight();
            List<int[]> endpoints = new ArrayList<>();
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    if (skelIP.get(x, y) != 0 && countNeighbors(skelIP, x, y, w, h) == 1) {
                        endpoints.add(new int[]{x, y});
                    }
                }
            }

            // Border fragments are only computed if a valid skeleton with two endpoints was found
            if (endpoints.size() >= 2) {
                int[] leftEndpoint = endpoints.get(0);
                for (int[] ep : endpoints) {
                    if (ep[0] < leftEndpoint[0]) leftEndpoint = ep;
                }

                FragmentData fragData = computeFragmentData(maskIP, pixelSizeUm, borderWidthUm, nFragments, leftEndpoint);
                List<List<int[]>> fragments = fragData.borderFragments;

                double[] gfpBorderSums = sumFragmentIntensities(fragments, gfpIP, gfpBackground);
                double[] dapiBorderSums = sumFragmentIntensities(fragments, dapiIP, dapiBackground);

                double[] gfpInteriorSums = sumInteriorIntensities(fragData.erodedMask, fragData.centroidX, fragData.centroidY,
                        fragData.startAngle, fragData.angleBoundaries, gfpIP, gfpBackground);
                double[] dapiInteriorSums = sumInteriorIntensities(fragData.erodedMask, fragData.centroidX, fragData.centroidY,
                        fragData.startAngle, fragData.angleBoundaries, dapiIP, dapiBackground);

                // Draw fragments with a distinct gray level per fragment index (1..nFragments)
                ImageProcessor fragIP = fragResult.getProcessor();
                for (int f = 0; f < fragments.size(); f++) {
                    int grayValue = f + 1;
                    for (int[] p : fragments.get(f)) {
                        fragIP.set(p[0], p[1], grayValue);
                    }
                }

                // Write one row per bacterium: DAPI/GFP total-intensity ratio for each fragment
                StringBuilder row = new StringBuilder(imgName + "\t" + time + "\t" + i);
                for (int f = 0; f < nFragments; f++) {
                    double borderRatio = (gfpBorderSums[f] != 0) ? dapiBorderSums[f] / gfpBorderSums[f] : Double.NaN;
                    double interiorRatio = (gfpInteriorSums[f] != 0) ? dapiInteriorSums[f] / gfpInteriorSums[f] : Double.NaN;
                    row.append("\t").append(dapiBorderSums[f]).append("\t").append(gfpBorderSums[f]).append("\t").append(borderRatio);
                    row.append("\t").append(dapiInteriorSums[f]).append("\t").append(gfpInteriorSums[f]).append("\t").append(interiorRatio);
                }
                row.append("\n");
                fragFile.write(row.toString());
                fragFile.flush();
            }

            ImageProcessor skelResultIP = skelResult.getProcessor();
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    if (skelIP.get(x, y) != 0) {
                        skelResultIP.set(x, y, 255);
                    }
                }
            }

            bact.changes = false;
            bact.close();
        }

        // Save the fragment label image separately, with a distinct-color LUT
        IJ.run(fragResult, "glasbey on dark", "");
        fragResult.setCalibration(img1.getCalibration());
        FileSaver fragFileSaver = new FileSaver(fragResult);
        fragFileSaver.saveAsTiff(outDir + imgName +"_fragments.tif");
        flush_close(fragResult);

        return skelResult;
    }


    /**
     * Compute the total (summed) pixel intensity of each fragment in a single channel image.
     *
     * @param fragments list of pixel-coordinate lists, one per fragment (from getBorderFragments)
     * @param channelIP channel image to measure
     * @return array of summed intensities, one value per fragment
     */

    public double[] sumFragmentIntensities(List<List<int[]>> fragments, ImageProcessor channelIP, double background) {
        int nFragments = fragments.size();
        double[] sums = new double[nFragments];

        for (int f = 0; f < nFragments; f++) {
            double sum = 0;
            for (int[] p : fragments.get(f)) {
                double corrected = channelIP.getPixelValue(p[0], p[1]) - background;
                if (corrected > 0) sum += corrected;
            }
            sums[f] = sum;
        }
        return sums;
    }

    /**
     * Sum background-corrected intensity over the interior of each pie-slice — i.e. the pixels
     * of the eroded mask (excluding the border ring, so no pixel is counted in both this and
     * sumFragmentIntensities) — using the same angular boundaries as the border fragments.
     */
    public double[] sumInteriorIntensities(ImageProcessor erodedMask, double centroidX, double centroidY,
                                           double startAngle, double[] angleBoundaries,
                                           ImageProcessor channelIP, double background) {

        int nFragments = angleBoundaries.length - 1;
        double[] sums = new double[nFragments];
        int w = erodedMask.getWidth(), h = erodedMask.getHeight();

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (erodedMask.get(x, y) == 0) continue; // not part of the interior

                double angle = Math.atan2(y - centroidY, x - centroidX);
                double relative = angle - startAngle;
                relative = ((relative % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);

                int fragIdx = nFragments - 1;
                for (int f = 0; f < nFragments; f++) {
                    if (relative >= angleBoundaries[f] && relative < angleBoundaries[f + 1]) {
                        fragIdx = f;
                        break;
                    }
                }

                double corrected = channelIP.getPixelValue(x, y) - background;
                if (corrected > 0) sums[fragIdx] += corrected;
            }
        }
        return sums;
    }

    private double estimateBackground(ImageProcessor channelIP, ImageProcessor labelIP) {
        int w = channelIP.getWidth();
        int h = channelIP.getHeight();

        List<Float> backgroundValues = new ArrayList<>();
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (labelIP.get(x, y) == 0) {
                    backgroundValues.add(channelIP.getPixelValue(x, y));
                }
            }
        }

        if (backgroundValues.isEmpty()) return 0;

        Collections.sort(backgroundValues);
        int mid = backgroundValues.size() / 2;
        if (backgroundValues.size() % 2 == 0) {
            return (backgroundValues.get(mid - 1) + backgroundValues.get(mid)) / 2.0;
        } else {
            return backgroundValues.get(mid);
        }
    }
    /**
     * Extend each skeleton endpoint in a straight line until it reaches the true edge
     * of the bacterium mask. The extension direction is estimated from the last few
     * skeleton pixels (PCA regression), which is more robust to local noise than using
     * only the endpoint's single neighbor.
     */
    private void extendSkeletonEndpoints(ImageProcessor skelIP, ImageProcessor maskIP) {
        int w = skelIP.getWidth();
        int h = skelIP.getHeight();
        int traceLength = 8; // number of skeleton pixels used to estimate the local direction

        List<int[]> endpoints = new ArrayList<>();
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (skelIP.get(x, y) == 0) continue;
                if (countNeighbors(skelIP, x, y, w, h) == 1) {
                    endpoints.add(new int[]{x, y});
                }
            }
        }

        for (int[] ep : endpoints) {
            List<int[]> path = traceBranch(skelIP, ep[0], ep[1], traceLength, w, h);
            if (path.size() < 2) continue;

            double[] dir = computeDirectionPCA(path, ep[0], ep[1]);
            if (dir == null) continue;

            double dx = dir[0];
            double dy = dir[1];
            double cx = ep[0];
            double cy = ep[1];

            while (true) {
                cx += dx;
                cy += dy;
                int px = (int) Math.round(cx);
                int py = (int) Math.round(cy);

                if (px < 0 || px >= w || py < 0 || py >= h) break;
                if (maskIP.get(px, py) == 0) break;

                skelIP.set(px, py, 255);
            }
        }
    }

    /**
     * Count the 8-connected non-zero neighbors of a pixel (used to detect skeleton endpoints,
     * which have exactly one neighbor).
     */
    private int countNeighbors(ImageProcessor ip, int x, int y, int w, int h) {
        int nCount = 0;
        for (int dyy = -1; dyy <= 1; dyy++) {
            for (int dxx = -1; dxx <= 1; dxx++) {
                if (dxx == 0 && dyy == 0) continue;
                int xx = x + dxx, yy = y + dyy;
                if (xx >= 0 && xx < w && yy >= 0 && yy < h && ip.get(xx, yy) != 0) {
                    nCount++;
                }
            }
        }
        return nCount;
    }


    /**
     * Walk the skeleton chain from an endpoint, pixel by pixel, up to maxSteps or until
     * a branch point (more than one unvisited neighbor) is reached. Returns the visited
     * path, starting at the endpoint.
     */
    private List<int[]> traceBranch(ImageProcessor skelIP, int startX, int startY, int maxSteps, int w, int h) {
        List<int[]> path = new ArrayList<>();
        boolean[][] visited = new boolean[w][h];

        int cx = startX, cy = startY;
        path.add(new int[]{cx, cy});
        visited[cx][cy] = true;

        for (int step = 0; step < maxSteps; step++) {
            int nextX = -1, nextY = -1;
            int candidates = 0;

            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int xx = cx + dx, yy = cy + dy;
                    if (xx >= 0 && xx < w && yy >= 0 && yy < h
                            && skelIP.get(xx, yy) != 0 && !visited[xx][yy]) {
                        candidates++;
                        nextX = xx;
                        nextY = yy;
                    }
                }
            }

            // Stop at a branch point rather than mixing directions from different branches
            if (candidates != 1) break;

            cx = nextX;
            cy = nextY;
            visited[cx][cy] = true;
            path.add(new int[]{cx, cy});
        }

        return path;
    }


    /**
     * Compute the principal direction of a skeleton segment via PCA (eigenvector of the
     * largest eigenvalue of the 2x2 covariance matrix), oriented to point from the segment
     * towards the given endpoint (outward direction).
     */
    private double[] computeDirectionPCA(List<int[]> path, int endX, int endY) {
        int n = path.size();
        if (n < 2) return null;

        double meanX = 0, meanY = 0;
        for (int[] p : path) {
            meanX += p[0];
            meanY += p[1];
        }
        meanX /= n;
        meanY /= n;

        double sxx = 0, syy = 0, sxy = 0;
        for (int[] p : path) {
            double dx = p[0] - meanX;
            double dy = p[1] - meanY;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }

        double theta = 0.5 * Math.atan2(2 * sxy, sxx - syy);
        double dx = Math.cos(theta);
        double dy = Math.sin(theta);

        double toEndX = endX - meanX;
        double toEndY = endY - meanY;
        if (dx * toEndX + dy * toEndY < 0) {
            dx = -dx;
            dy = -dy;
        }

        double norm = Math.sqrt(dx * dx + dy * dy);
        return new double[]{dx / norm, dy / norm};
    }


    /**
     * Manual binary erosion (independent of ImageJ's global "black background" setting).
     * iterations = number of pixel layers removed from the mask border.
     */
    private ByteProcessor erode(ImageProcessor maskIP, int iterations) {
        ByteProcessor current = (ByteProcessor) maskIP.convertToByte(false).duplicate();
        int w = current.getWidth();
        int h = current.getHeight();

        for (int it = 0; it < iterations; it++) {
            ByteProcessor next = (ByteProcessor) current.duplicate();
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    if (current.get(x, y) == 0) continue;
                    boolean erodePixel = false;
                    for (int dy = -1; dy <= 1 && !erodePixel; dy++) {
                        for (int dx = -1; dx <= 1 && !erodePixel; dx++) {
                            if (dx == 0 && dy == 0) continue;
                            int xx = x + dx, yy = y + dy;
                            if (xx < 0 || xx >= w || yy < 0 || yy >= h || current.get(xx, yy) == 0) {
                                erodePixel = true;
                            }
                        }
                    }
                    if (erodePixel) next.set(x, y, 0);
                }
            }
            current = next;
        }
        return current;
    }


    /**
     * Using the centroid angle avoids the ambiguity issues of contour- or skeleton-based
     * ordering (which are loops/paths where two spatially close pixels can have very
     * different positions along the path near thin tips). The angle around a single fixed
     * point is a strictly monotonic function of position, so it never causes fragment mixing.
     *
     * @param maskIP        original (pre-erosion) binary mask of the bacterium
     * @param pixelSizeUm   pixel size in µm (assumes isotropic calibration)
     * @param borderWidthUm desired border thickness in µm
     * @param nFragments    number of fragments requested
     * @param leftEndpoint  skeleton's left endpoint {x, y}, used to fix the sector origin angle
     */
    /**
     * Holds the border fragments together with the angular boundaries and eroded mask used
     * to build them, so the same boundaries can be reused to sum intensities over each
     * fragment's interior (same angular sector, but inside the eroded mask instead of the
     * border ring).
     */
    private static class FragmentData {
        List<List<int[]>> borderFragments;
        double[] angleBoundaries; // size nFragments+1, relative angles in [0, 2*PI], increasing
        double centroidX, centroidY;
        double startAngle;
        ByteProcessor erodedMask;
    }

    /**
     * Divide the bacterium border (original mask minus eroded mask) into nFragments of equal
     * pixel count, as angular sectors ("pie slices") around the bacterium centroid. The
     * sector boundary starts at the angle of the skeleton's left endpoint (its "12 o'clock").
     */
    public List<List<int[]>> getBorderFragments(ImageProcessor maskIP, double pixelSizeUm,
                                                double borderWidthUm, int nFragments, int[] leftEndpoint) {
        return computeFragmentData(maskIP, pixelSizeUm, borderWidthUm, nFragments, leftEndpoint).borderFragments;
    }

    /**
     * Compute the border fragments and the angular sector boundaries used to build them.
     * Reusing these boundaries lets the interior of each sector be measured separately
     * (see sumInteriorIntensities) without recomputing the geometry.
     */
    private FragmentData computeFragmentData(ImageProcessor maskIP, double pixelSizeUm,
                                             double borderWidthUm, int nFragments, int[] leftEndpoint) {

        int erosionPixels = Math.max(1, (int) Math.round(borderWidthUm / pixelSizeUm));
        ByteProcessor eroded = erode(maskIP, erosionPixels);

        int w = maskIP.getWidth();
        int h = maskIP.getHeight();

        FragmentData data = new FragmentData();
        data.borderFragments = new ArrayList<>();
        for (int f = 0; f < nFragments; f++) data.borderFragments.add(new ArrayList<>());
        data.angleBoundaries = new double[nFragments + 1];
        data.erodedMask = eroded;

        // Centroid of the ORIGINAL mask (before erosion) is the reference center for angles
        double cx = 0, cy = 0;
        int maskCount = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (maskIP.get(x, y) != 0) {
                    cx += x;
                    cy += y;
                    maskCount++;
                }
            }
        }
        if (maskCount == 0) return data;
        cx /= maskCount;
        cy /= maskCount;
        data.centroidX = cx;
        data.centroidY = cy;

        // Reference angle = angle of the skeleton's left endpoint, seen from the centroid
        double startAngle = Math.atan2(leftEndpoint[1] - cy, leftEndpoint[0] - cx);
        data.startAngle = startAngle;

        List<double[]> borderPositions = new ArrayList<>(); // {x, y, relativeAngle}
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (maskIP.get(x, y) == 0 || eroded.get(x, y) != 0) continue; // not a border pixel

                double angle = Math.atan2(y - cy, x - cx);
                double relative = angle - startAngle;
                relative = ((relative % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);

                borderPositions.add(new double[]{x, y, relative});
            }
        }

        int total = borderPositions.size();
        if (total < nFragments) {
            System.out.println("Warning: only " + total + " border pixels for " + nFragments + " requested fragments.");
        }

        borderPositions.sort((a, b) -> Double.compare(a[2], b[2]));

        int baseSize = total / nFragments;
        int remainder = total % nFragments;
        int idx = 0;

        data.angleBoundaries[0] = 0;
        data.angleBoundaries[nFragments] = 2 * Math.PI;

        for (int f = 0; f < nFragments; f++) {
            if (f > 0) {
                data.angleBoundaries[f] = (idx < total) ? borderPositions.get(idx)[2] : data.angleBoundaries[f - 1];
            }
            int fragSize = baseSize + (f < remainder ? 1 : 0);
            for (int k = 0; k < fragSize && idx < total; k++, idx++) {
                double[] p = borderPositions.get(idx);
                data.borderFragments.get(f).add(new int[]{(int) p[0], (int) p[1]});
            }
        }

        return data;
    }


    /**
     * Compute and write bacteria/DNA parameters (surface, length, DNA count, DNA intensity,
     * distance to bacterium center, volume ratio) to the results file.
     */
    public void saveResults(Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop,
                            ImagePlus dnaImg, String imgName, int time, BufferedWriter file) throws IOException {

        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            int bactLabel = (int) bact.getLabel();
            double bactSurf = new MeasureVolume(bact).getVolumeUnit();
            double bactLength = new MeasureFeret(bact).getFeret1Unit()
                    .distance(new MeasureFeret(bact).getFeret2Unit()) * cal.pixelWidth;

            // ADN liés à 1 seule bactérie
            List<Object3DInt> dnaBactList = dnaPop.getObjects3DInt().stream()
                    .filter(dna -> (int) dna.getIdObject() == bactLabel)
                    .collect(Collectors.toList());

            // ADN liés à plusieurs bactéries dont celle-ci
            Objects3DIntPopulation dnaBactPop = findDnaInBact(bactLabel, dnaPop);
            int dnaNb = dnaBactPop.getNbObjects();

            if (dnaNb == 0) {
                file.write(imgName + "\t" + time + "\t" + bactLabel + "\t" + bactSurf + "\t" + bactLength + "\t" + dnaNb + "\n");
            } else {
                for (Object3DInt dna : dnaBactPop.getObjects3DInt()) {
                    double dnaSurf = new MeasureVolume(dna).getVolumeUnit();
                    double dnaInt = new MeasureIntensity(dna, ImageHandler.wrap(dnaImg))
                            .getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                    double dnaDist = dnaBactDistance(dna, bact);
                    double totalVol = dnaTotalBactVol.getOrDefault((int) dna.getIdObject(), bactSurf);

                    file.write(imgName + "\t" + time + "\t" + bactLabel + "\t" + bactSurf + "\t" + bactLength + "\t"
                            + dnaNb + "\t" + dna.getLabel() + "\t" + dnaSurf + "\t" + dnaInt + "\t" + dnaDist + "\t" + totalVol + "\t" + (dnaSurf / totalVol) + "\t" + "\n");
                }
            }
            file.flush();
        }
    }


    /**
     * Compute the distance between a DNA object's centroid and its parent bacterium's centroid.
     */
    private double dnaBactDistance(Object3DInt dna, Object3DInt bact) {
        Voxel3D bactCenter = new MeasureCentroid(bact).getCentroidAsVoxel();
        Voxel3D dnaCenter = new MeasureCentroid(dna).getCentroidAsVoxel();
        return bactCenter.distance(dnaCenter) * cal.pixelWidth;
    }


    /**
     * Save labeled bacteria and DNA overlays merged with their respective intensity channel,
     * as two TIFF files.
     */
    public void drawResults(ImagePlus img1, ImagePlus img2, Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop, String imgName, String outDir) {
        ImageHandler imgBact = ImageHandler.wrap(img1).createSameDimensions();
        bactPop.drawInImage(imgBact);
        IJ.run(imgBact.getImagePlus(), "glasbey on dark", "");
        ImagePlus[] imgColors1 = {imgBact.getImagePlus(), null, null, img1};
        ImagePlus imgOut1 = new RGBStackMerge().mergeHyperstacks(imgColors1, false);
        imgOut1.setCalibration(cal);
        new FileSaver(imgOut1).saveAsTiff(outDir + imgName + "_bacteria.tif");

        ImageHandler imgDna = ImageHandler.wrap(img2).createSameDimensions();
        dnaPop.drawInImage(imgDna);
        ImagePlus[] imgColors2 = {imgDna.getImagePlus(), null, null, img2};
        ImagePlus imgOut2 = new RGBStackMerge().mergeHyperstacks(imgColors2, false);
        imgOut2.setCalibration(cal);
        new FileSaver(imgOut2).saveAsTiff(outDir + imgName + "_DNA.tif");

        flush_close(imgBact.getImagePlus());
        flush_close(imgDna.getImagePlus());
        flush_close(imgOut1);
        flush_close(imgOut2);
    }
}