import GFP_density_border_Tools.Tools;
import ij.*;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.scijava.util.ArrayUtils;


/**
 * Detect bacteria and DNA in them with Omnipose
 * @author Orion-CIRB
 */
public class GFP_density_border implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public String outDirResults = "";
    public BufferedWriter results;
   
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            } 
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }  
            
            // Find images with extension
            String file_ext = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, file_ext);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
             String header = "Image name\tTime\t# bacterium\tBacterium surface (µm2)\tBacterium length (µm)\tDNA number\t# DNA\tDNA surface (µm2)\tDNA total intensity\t"
                     + "DNA center to bacterium center (µm)\tTotal Bacterium Surface\t Ratio \n";
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();

            // Avant la boucle sur les images (là où "results" est probablement déjà ouvert)
            BufferedWriter fragResults = new BufferedWriter(new FileWriter(outDirResults + "fragmentResults.xls"));
            tools.writeFragmentHeader(fragResults, tools.nFragments); // 10 = nFragments, à garder cohérent avec l'appel plus bas
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channels name
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Dialog box
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showMessage("Error", "Plugin canceled");
                return;
            }
            
            for (String f : imageFiles) {
                reader.setId(f);              
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setSplitChannels(true);
                
                int series = reader.getSeriesCount();
                for (int s = 0; s < series; s++) {
                    reader.setSeries(s);
                    options.setSeriesOn(s, true);
                    String seriesName = meta.getImageName(s);
                    
                    int time = reader.getSizeT();
                    for (int t = 0; t < time; t++) {
                        tools.print("--- ANALYZING IMAGE " + seriesName + " at time " + (t+1) + " ---");
                        options.setTBegin(s, t);
                        options.setTEnd(s, t);

                        // Open bacteria channel
                        int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                        System.out.println("Opening phase channel "+chs[0] );
                        ImagePlus bactStack = BF.openImagePlus(options)[indexCh];
                        ImagePlus imgBact = tools.doZProjection(bactStack, ZProjector.AVG_METHOD);
                        tools.flush_close(bactStack);

                        // Open Gfp chan
                        indexCh = ArrayUtils.indexOf(channels, chs[2]);
                        ImagePlus GFPStack = BF.openImagePlus(options)[indexCh];
                        ImagePlus imgGFP = tools.doZProjection(GFPStack, ZProjector.AVG_METHOD);
                        tools.flush_close(GFPStack);

                        // Detect bacteria with Omnipose
                        tools.print("- Detecting bacteria -");
                        Objects3DIntPopulation bactPop = tools.omniposeDetection(imgBact, tools.omniposeBactModel, tools.minBactSurface, tools.maxBactSurface, true);
                        System.out.println(bactPop.getNbObjects() + " bacteria found");

                        // Open DNA channel
                        indexCh = ArrayUtils.indexOf(channels, chs[1]);
                        System.out.println("Opening DNA channel "+chs[1]);
                        ImagePlus dnaStack = BF.openImagePlus(options)[indexCh];
                        ImagePlus imgDna = tools.doZProjection(dnaStack, ZProjector.AVG_METHOD);
                        tools.flush_close(dnaStack);

                        // Detect DNA with Omnipose
                        tools.print("- Detecting DNA -");
                        Objects3DIntPopulation dnaPop = tools.omniposeDetection(imgDna, tools.omniposeDnaModel, tools.minDnaSurface, tools.maxDnaSurface, false);
                        System.out.println(dnaPop.getNbObjects() + " DNA found");
                        tools.dnaBactLink(bactPop, dnaPop);
                        tools.dnaBactNumber(bactPop,dnaPop);


                        System.out.println(dnaPop.getNbObjects() + " DNA found in bacteria");

                        // Save results
                        tools.print("- Saving results -");
                        tools.saveResults(bactPop, dnaPop, imgDna, seriesName, t+1, results);

                        // Save images
                        ImagePlus img_skel = tools.bactBorderIntensities(imgBact, bactPop, imgGFP, imgDna, tools.borderWidth, tools.nFragments,
                                seriesName + "_t" + (t+1), t+1, fragResults, outDirResults);
                        tools.drawResults(img_skel, imgDna, bactPop, dnaPop, seriesName+"_t"+(t+1), outDirResults);

                        tools.flush_close(imgBact);
                        tools.flush_close(imgDna);
                    }
                    options.clearSeries();
                    tools.initImageVariables();
                }
            }

            fragResults.close();
            tools.print("--- All done! ---");
            
        }   catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(GFP_density_border.class.getName()).log(Level.SEVERE, null, ex);
        }  
    }
}    
