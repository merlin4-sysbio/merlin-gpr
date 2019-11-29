/**
 * 
 */
package pt.uminho.ceb.biosystems.merlin.gpr.rules.operations;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.core.sequence.template.AbstractSequence;

import es.uvigo.ei.aibench.core.operation.annotation.Cancel;
import es.uvigo.ei.aibench.core.operation.annotation.Direction;
import es.uvigo.ei.aibench.core.operation.annotation.Operation;
import es.uvigo.ei.aibench.core.operation.annotation.Port;
import es.uvigo.ei.aibench.core.operation.annotation.Progress;
import es.uvigo.ei.aibench.workbench.Workbench;
import pt.uminho.ceb.biosystems.merlin.aibench.datatypes.WorkspaceAIB;
import pt.uminho.ceb.biosystems.merlin.aibench.gui.CustomGUI;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.LoadFromConf;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.MerlinUtils;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.TimeLeftProgress;
import pt.uminho.ceb.biosystems.merlin.core.containers.gpr.ReactionsGPR_CI;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.Method;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SequenceType;
import pt.uminho.ceb.biosystems.merlin.gpr.rules.core.FilterModelReactions;
import pt.uminho.ceb.biosystems.merlin.gpr.rules.core.IdentifyGenomeSubunits;
import pt.uminho.ceb.biosystems.merlin.local.alignments.core.AlignmentsUtils;
import pt.uminho.ceb.biosystems.merlin.processes.WorkspaceProcesses;
import pt.uminho.ceb.biosystems.merlin.services.ProjectServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelGenesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelSequenceServices;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;

/**
 * @author ODias
 *
 */
@Operation(description="Generate and integrate gene-reactions connections.",name="GPRs generator.")
public class CreateGPRs implements PropertyChangeListener {

	private long reference_organism_id;
	private double similarity_threshold;
	private double referenceTaxonomyThreshold;
	private boolean compareToFullGenome;
	private boolean identifyGPRs;
	private boolean integrateToDatabase;
	private boolean keepReactionsWithNotes;
	private boolean isCompartimentalised;
	private boolean generateGPRs;
	private TimeLeftProgress  progress = new TimeLeftProgress();
	private Map<String, AbstractSequence<?>> genome;
	private boolean keepManualReactions;
	private boolean removeReactions;
	private AtomicBoolean cancel = new AtomicBoolean(false);
	private Method method;
	private long startTime;
	private int dataSize;
	private WorkspaceAIB project;
	private float threshold;

	/**
	 * @param project
	 */
	@Port(name="workspace",description="select workspace",direction=Direction.INPUT,order=1, validateMethod="checkProject")

	public void setProject(WorkspaceAIB project) {
		
		this.startTime = GregorianCalendar.getInstance().getTimeInMillis();
		
		this.reference_organism_id = project.getTaxonomyID();
		try {
			this.isCompartimentalised = ProjectServices.isCompartmentalisedModel(this.project.getName());
		}
		catch (Exception e1) {
			
			Workbench.getInstance().error(e1);
		}

		try {

			Map<String, String> settings = LoadFromConf.loadGPRsettings(FileUtils.getConfFolderPath());

			this.similarity_threshold = Float.parseFloat(settings.get("similarity_threshold"));
			this.referenceTaxonomyThreshold = Float.parseFloat(settings.get("referenceTaxonomyThreshold"));
			this.compareToFullGenome = Boolean.parseBoolean(settings.get("compareToFullGenome"));
			this.identifyGPRs = Boolean.parseBoolean(settings.get("identifyGPRs"));
			this.generateGPRs = Boolean.parseBoolean(settings.get("generateGPRs"));
			this.keepReactionsWithNotes = Boolean.parseBoolean(settings.get("keepReactionsWithNotes"));
			this.keepManualReactions = Boolean.parseBoolean(settings.get("keepManualReactions"));
			this.integrateToDatabase = Boolean.parseBoolean(settings.get("integrateToDatabase"));
			this.threshold = Float.parseFloat(settings.get("threshold"));
			this.removeReactions = Boolean.parseBoolean(settings.get("removeReactions"));

			if(settings.get("method").equalsIgnoreCase("blast"))
				this.method = Method.Blast;
			else
				this.method = Method.SmithWaterman;

			
			if(!this.method.equals(Method.Blast) || AlignmentsUtils.checkBlastInstalation()){

				boolean identifiedWithoutErros = false;
				if(!this.cancel.get()) {
					
					if(this.identifyGPRs && !this.cancel.get()) {

						Map<String, List<String>> ec_numbers = ModelGenesServices.getGPRsECNumbers(this.project.getName(), this.isCompartimentalised);

						IdentifyGenomeSubunits identifyGenomeSubunits = new IdentifyGenomeSubunits(ec_numbers, genome, reference_organism_id, similarity_threshold, 
								referenceTaxonomyThreshold, this.method, compareToFullGenome);

						String wsTaxonomyGprsFolderPath = FileUtils.getWorkspaceTaxonomyGprsFolderPath(project.getName(), project.getTaxonomyID());
						identifyGenomeSubunits.setWsTaxonomyFolderPath(wsTaxonomyGprsFolderPath);
						identifyGenomeSubunits.addPropertyChangeListener(this);
						identifyGenomeSubunits.setCancel(this.cancel);
						identifiedWithoutErros = identifyGenomeSubunits.runIdentification(false, project.getName());

					}
					
				}
				
				if(!this.cancel.get()) {
					
					if(!this.identifyGPRs || identifiedWithoutErros) {

						if(this.generateGPRs && !this.cancel.get()) {

							Map<String, ReactionsGPR_CI> ret = IdentifyGenomeSubunits.runGPRsAssignment(this.project.getName(),this.threshold);

							FilterModelReactions f = new FilterModelReactions(this.project.getName(), this.isCompartimentalised);
							f.filterReactions(ret);

							if(this.integrateToDatabase && !this.cancel.get()) {

								if(this.removeReactions)
									f.removeReactionsFromModel(keepReactionsWithNotes, this.keepManualReactions);
								
								f.setModelGPRsFromTool();
							}
						}

						if(this.cancel.get()) {

							Workbench.getInstance().warn("gene-protein-reactions rules finding canceled");
						}
						else {

							MerlinUtils.updateAllViews(project.getName());
							Workbench.getInstance().info("gene-protein-reactions rules finding finished");
						}
					}
				}
				else {

					if(this.cancel.get())
						Workbench.getInstance().warn("gene-protein-reactions rules finding canceled");

					else if(!this.identifyGPRs)
						Workbench.getInstance().error("merlin found some errors whilst performing this operation, please try again later");
				}
			}
			else{
				Workbench.getInstance().error("blast+ not found! Please install it or change 'method' parameter in\n"
						+ "\"\\conf\\gpr_settings.conf\" file to 'smithwaterman' to run gprs identification.");
			}
		} 
		catch (Exception e) {

			Workbench.getInstance().error("error "+e.getMessage()+" has occured.");
			e.printStackTrace();
		}

	}

	/**
	 * @param similarity_threshold
	 */
	public void checkDouble(double double_) {

		if(double_>1 || double_<0) {

			throw new IllegalArgumentException("please set a valid double (0<double<1).");
		}
	}

	/**
	 * @return the progress
	 */
	@Progress
	public TimeLeftProgress getProgress() {

		return progress;
	}

	/**
	 * @param cancel the cancel to set
	 */
	@Cancel
	public void setCancel() {
		
		String[] options = new String[2];
		options[0] = "yes";
		options[1] = "no";
		
		int result = CustomGUI.stopQuestion("cancel confirmation", "are you sure you want to cancel the operation?", options);
		
		if (result == 0) {
			
			progress.setTime(0, 0, 0);
			this.cancel.set(true);
			
		}
	}

	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		
		if(evt.getPropertyName().equalsIgnoreCase("size"))
			this.dataSize = new AtomicInteger((int) evt.getNewValue()).get();
		
		if(evt.getPropertyName().equalsIgnoreCase("sequencesCounter")) {
			int sequencesCounter = (int) evt.getNewValue();
			this.progress.setTime((GregorianCalendar.getInstance().getTimeInMillis() - startTime), sequencesCounter, dataSize);
		}
	}
	
	/**
	 * @param project
	 */
	public void checkProject(WorkspaceAIB project) {

		if(project == null) {

			throw new IllegalArgumentException("No ProjectGUISelected!");
		}
		else {
			//			String dbName = project.getName();
			//			Long taxID = project.getTaxonomyID();
			try {
				
				this.project = project;

				if(!ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.PROTEIN)){//!Project.isFaaFiles(dbName,taxID) && !Project.isFnaFiles(dbName,taxID)) {

					throw new IllegalArgumentException("please set the genome fasta file ('.faa') to perform gene-protein-reactions rules identification");
				}
				else if(project.getTaxonomyID()<0) {

					throw new IllegalArgumentException("please enter the organism taxonomic identification from NCBI taxonomy to perform this operation");
				}
				else {

					try {

						this.genome = ModelSequenceServices.getGenomeFromDatabase(this.project.getName(), SequenceType.PROTEIN);

						if(this.identifyGPRs && (this.genome==null || this.genome.isEmpty())) {
							WorkspaceProcesses.createFaaFile(this.project.getName(), this.project.getTaxonomyID()); // method creates ".faa" files only if they do not exist
						}
					} 
					catch (Exception e) {

						e.printStackTrace();
						throw new IllegalArgumentException("please set the project fasta files");
					}
				}

			} catch(Exception e1) {

				e1.printStackTrace();
			}
		}
	}


}
