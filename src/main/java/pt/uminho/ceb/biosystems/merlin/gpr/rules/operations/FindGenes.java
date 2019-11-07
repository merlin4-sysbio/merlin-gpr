package pt.uminho.ceb.biosystems.merlin.gpr.rules.operations;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.nbio.core.sequence.template.AbstractSequence;

import es.uvigo.ei.aibench.core.operation.annotation.Direction;
import es.uvigo.ei.aibench.core.operation.annotation.Operation;
import es.uvigo.ei.aibench.core.operation.annotation.Port;
import es.uvigo.ei.aibench.core.operation.annotation.Progress;
import es.uvigo.ei.aibench.workbench.Workbench;
import pt.uminho.ceb.biosystems.merlin.aibench.datatypes.model.ModelReactionsAIB;
import pt.uminho.ceb.biosystems.merlin.aibench.gui.CustomGUI;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.LoadFromConf;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.TimeLeftProgress;
import pt.uminho.ceb.biosystems.merlin.core.containers.alignment.AlignmentContainer;
import pt.uminho.ceb.biosystems.merlin.core.containers.model.GeneContainer;
import pt.uminho.ceb.biosystems.merlin.core.datatypes.WorkspaceDataTable;
import pt.uminho.ceb.biosystems.merlin.core.datatypes.WorkspaceGenericDataTable;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.Method;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SequenceType;
import pt.uminho.ceb.biosystems.merlin.gpr.rules.core.IdentifyGenomeSubunits;
import pt.uminho.ceb.biosystems.merlin.gpr.rules.gui.FillGapReaction;
import pt.uminho.ceb.biosystems.merlin.local.alignments.core.AlignmentsUtils;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelEnzymesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelGenesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelSequenceServices;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

@Operation(name="FindGenes", description="find genes for gap reactions")
public class FindGenes implements PropertyChangeListener {

	private Integer reactionID;
	private Map<String, List<String>> ecNumbers;
	private int dataSize;
	private TimeLeftProgress  progress = new TimeLeftProgress();
	private long startTime;

	@Port(direction=Direction.INPUT, name="reactionID", description = "", order=1)
	public void setReactionID(Integer reactionID){

		this.reactionID = reactionID;
	};

	@Port(direction=Direction.INPUT, name="modelReactions", description = "" , order=2)
	public void setReactionsInterface(ModelReactionsAIB reactionInterface) throws Exception{

 		this.startTime = GregorianCalendar.getInstance().getTimeInMillis();

		long taxonomyID = reactionInterface.getWorkspace().getTaxonomyID();
		String databaseName = reactionInterface.getWorkspace().getName();

		Map<String, String> credentials = LoadFromConf.loadReactionsThresholds(FileUtils.getConfFolderPath());

		double similarityThreshold = Double.valueOf(credentials.get("similarity_threshold"));
		Double reference_taxo_threshold = Double.valueOf(credentials.get("reference_taxo_threshold"));

		int flag = -1;
		WorkspaceDataTable data = null;

		Pair<Integer, WorkspaceDataTable> pair = FillGapReaction.openFile(databaseName, taxonomyID, this.reactionID, similarityThreshold);
		flag = pair.getA();
		data = pair.getB();
		
		this.ecNumbers = new HashMap<String, List<String>>();
		
		Method method = Method.Blast;
		boolean go = true;

		if(!AlignmentsUtils.checkBlastInstalation()){

			method = Method.SmithWaterman;
			go = continueIdentification();
		}

		if(go){

			if(flag == -1){

				List<String> ecs = ModelEnzymesServices.getEcNumbersList(reactionInterface.getWorkspace().getName(), reactionID);

				for(String ec : ecs)
					this.ecNumbers.put(ec, new ArrayList<>());
				 
				Map<String, AbstractSequence<?>> genome =  ModelSequenceServices.getGenomeFromDatabase(reactionInterface.getWorkspace().getName(), SequenceType.PROTEIN);

				IdentifyGenomeSubunits i = new IdentifyGenomeSubunits(ecNumbers, genome, taxonomyID, similarityThreshold, 
						reference_taxo_threshold, method , true);

//				String wsTaxonomyFolderPath = FileUtils.getWorkspaceTaxonomyFolderPath(reactionInterface.getWorkspace().getName(), reactionInterface.getWorkspace().getTaxonomyID());
				String wsTaxonomyGprsFolderPath = FileUtils.getWorkspaceTaxonomyFillGapsFolderPath(reactionInterface.getWorkspace().getName(), reactionInterface.getWorkspace().getTaxonomyID());

				i.setWsTaxonomyFolderPath(wsTaxonomyGprsFolderPath);
//				i.setWsTaxonomyTempFolderPath(wsTaxonomyGprsFolderPath);
				i.setCancel(new AtomicBoolean(false));
				boolean identifiedWithoutErros = i.runIdentification(true, databaseName);

				if(identifiedWithoutErros) {
					
					ConcurrentLinkedQueue<AlignmentContainer> alignmentResults=i.findGapsResult();
					data = this.createTable(databaseName, alignmentResults);
					this.writeNewFile(databaseName, taxonomyID, reactionID, similarityThreshold, data);
				}
				else{ 

					Workbench.getInstance().warn("an error ocurred");
				}
			}

			new FillGapReaction(reactionInterface, reactionID, similarityThreshold);
		}
		else {
			
			Workbench.getInstance().warn("find genes job canceled");
		}

	};



	/**
	 * @param databaseName
	 * @param alignmentResults
	 * @return
	 */
	public WorkspaceDataTable createTable(String databaseName, ConcurrentLinkedQueue<AlignmentContainer> alignmentResults) {
		
		List<String> alignCol = new ArrayList<>();
		alignCol.add("locus tag");
		alignCol.add("sequence id");
		alignCol.add("orthologous");
		alignCol.add("alignment score");
		alignCol.add("ec number");
		alignCol.add("query coverage");
		alignCol.add("target's coverage");

		WorkspaceDataTable data = new WorkspaceGenericDataTable(alignCol, "alignment table","alignment results");

		try {
			for (AlignmentContainer alignmentContainer : alignmentResults) {		

				ArrayList<String> line = new ArrayList<>();

				if(alignmentContainer.getQuery() != null || alignmentContainer.getAlignmentScore() != -1.0){

					GeneContainer container = ModelGenesServices.getGeneByQuery(databaseName, alignmentContainer.getTarget());
					if(container == null)
						line.add(alignmentContainer.getTarget());
					else
						line.add(container.getLocusTag());
					line.add(alignmentContainer.getTarget());
					line.add(alignmentContainer.getQuery());
					line.add(alignmentContainer.getAlignmentScore()+"");
					String ec = this.ecNumbers.keySet().toString();
					line.add(ec.replace("[", "").replace("]", "").trim());
					line.add(alignmentContainer.getCoverageQuery()+"");
					line.add(alignmentContainer.getCoverageTarget()+"");
					data.addLine(line);

				}
			}
		} catch (Exception e) {
			Workbench.getInstance().error(e);
			e.printStackTrace();
		}
		return data;

	}

	/**
	 * @param databaseName
	 * @param taxonomyID
	 * @param reactionID
	 * @param similarityThreshold
	 * @param data
	 * @throws IOException
	 */
	private void writeNewFile(String databaseName, long taxonomyID, int reactionID, double similarityThreshold, WorkspaceDataTable data) throws IOException{

		String path = FileUtils.getWorkspaceTaxonomyFillGapsFolderPath(databaseName, taxonomyID) + "FillGapReactions.txt";

		try {
			File file = new File(path);

			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(file,true)));

			if(data.getRowCount() > 0){

				for (int i=0; i < data.getRowCount(); i++){

					out.write(reactionID+"\t"+similarityThreshold+"\t");

					for (int j=0; j < data.getColumnCount(); j++){
						out.write(data.getValueAt(i,j)+"\t");
					}
					out.write("\n");
				}
			}
			else{
				out.write(reactionID+"\t"+similarityThreshold+"\t"+"empty"+"\n");
			}
			out.close();

		} 
		catch (Exception e) {
			e.printStackTrace();
		}

	}


	public boolean continueIdentification(){

		int i =CustomGUI.stopQuestion("Find genes alignment method", "Blast+ not found! SmithWaterman alignment method will be used instead. Do you wish to continue?",
				new String[]{"Continue", "Abort", "Info"});

		if(i<2) {

			switch (i)
			{
			case 0:return true;
			case 1: return false;
			default:return true;
			}
		}
		else{
			Workbench.getInstance().info("merlin has detected that the default alignment search tool, NCBI BLAST+, is not installed on your computer.\n"
					+ "For this reason 'find genes' operation will use SmithWaterman alignment method to perform similarity searches. Note that this\n"
					+ "method will be considerably slower then BLAST.\n"
					+ "In order to use BLAST, enter 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/', download and install NCBI BLAST+.\n"
					+ "If the problem persists after the instalation you may need to add BLAST instalation directory to windows class path.\n"
					+ "For further information please contact us to support@merlin-sysbio.org");

			return continueIdentification();
		}

	}

	@Progress(progressDialogTitle = "find genes", modal = false, workingLabel = "finding genes for reaction", preferredWidth = 300, preferredHeight=200)
	public TimeLeftProgress getProgress(){

		return progress;
	}


	@Override
	public void propertyChange(PropertyChangeEvent evt) {

		if(evt.getPropertyName().equalsIgnoreCase("size"))
			this.dataSize = (int) evt.getNewValue();

		if(evt.getPropertyName().equalsIgnoreCase("sequencesCounter")) {
			int sequencesCounter = (int) evt.getNewValue();
			this.progress.setTime((GregorianCalendar.getInstance().getTimeInMillis() - startTime), sequencesCounter, dataSize);
		}
	}

}
