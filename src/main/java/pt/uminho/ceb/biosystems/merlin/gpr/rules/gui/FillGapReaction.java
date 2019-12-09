package pt.uminho.ceb.biosystems.merlin.gpr.rules.gui;


import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

import es.uvigo.ei.aibench.workbench.Workbench;
import es.uvigo.ei.aibench.workbench.utilities.Utilities;
import pt.uminho.ceb.biosystems.merlin.gui.datatypes.model.ModelReactionsAIB;
import pt.uminho.ceb.biosystems.merlin.gui.utilities.CreateImageIcon;
import pt.uminho.ceb.biosystems.merlin.gui.utilities.LoadFromConf;
import pt.uminho.ceb.biosystems.merlin.gui.utilities.MerlinUtils;
import pt.uminho.ceb.biosystems.merlin.core.datatypes.WorkspaceDataTable;
import pt.uminho.ceb.biosystems.merlin.core.datatypes.WorkspaceGenericDataTable;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SourceType;
import pt.uminho.ceb.biosystems.merlin.services.ProjectServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelEnzymesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelGenesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelReactionsServices;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;

public class FillGapReaction extends JDialog {

	private static final long serialVersionUID = -1L;
	private long reference_organism_id;
	private ModelReactionsAIB reactionInterface;

	private WorkspaceDataTable data;
	private int reactionID;
	private JButton jButtonOK, jButtonCancel;

	/**
	 * @param reactionInterface
	 * @param reactionID
	 * @param similarityThreshold
	 * @param statement
	 */
	public FillGapReaction(ModelReactionsAIB reactionInterface, int reactionID, double similarityThreshold) {

		super(Workbench.getInstance().getMainFrame());
		this.reactionID=reactionID;
		this.reactionInterface = reactionInterface;
		Utilities.centerOnOwner(this);
		this.reference_organism_id = reactionInterface.getWorkspace().getTaxonomyID();

		String databaseName = reactionInterface.getWorkspace().getName();

		int flag = -1;

		Pair<Integer, WorkspaceDataTable> pair = openFile(databaseName, reference_organism_id, this.reactionID, similarityThreshold);
		flag = pair.getA();
		this.data = pair.getB();

		if(flag==1){
			initGui();
			Utilities.centerOnOwner(this);
			this.setAlwaysOnTop(true);
			this.toFront(); 
			this.setVisible(true);

		}
		else { 
			Workbench.getInstance().warn("nothing found");
		} 
	}

	/**
	 * @throws Exception 
	 * 
	 */
	private void apply(List<String> entries) throws Exception {


		Map<String, List<Integer>> reactions  = new HashMap<>();
		for(int i = 0; i < entries.size(); i++) {

			int geneID = ModelGenesServices.loadGene(this.reactionInterface.getWorkspace().getName(), null, entries.get(i), null, null, null, null, SourceType.KO);

			String ecs = ((String) data.getValueAt(i, 4));

			String[] ecNumbers = ecs.split(",");

			Set<String> ec = new HashSet<String>();


			for(String ecNumber : ecNumbers) {
				ec.add(ecNumber.trim());				
			}
			reactions.putAll(ModelEnzymesServices.loadEnzymeGetReactions(this.reactionInterface.getWorkspace().getName(),geneID, ec, null, 
					true, true, false, ProjectServices.isCompartmentalisedModel(reactionInterface.getWorkspace().getName())));	

		}

		if (reactions.keySet().isEmpty()){
			Workbench.getInstance().warn("no reactions to be included in model");
		}
		else {
			Set<String> react  = new HashSet<String>();

			for(String key : reactions.keySet()){

				List<Integer> identifiers = reactions.get(key);

				try {

					react.addAll(ModelReactionsServices.getReactionsNames(reactionInterface.getWorkspace().getName(), identifiers).values());

				} catch (Exception e) {
					e.printStackTrace();
				}

			}

			windowToBack();

			if(react.size() > 0){

				String reacts = react.toString();
				Workbench.getInstance().info("the following reactions were added to the model:\n"+reacts.replace("[", "").replace("]", ""));
				String workspace = reactionInterface.getWorkspace().getName();
				MerlinUtils.updateMetabolicViews(workspace);

			}
			else	
				Workbench.getInstance().warn("no reactions to be included in model");

		}
	}

	public static Pair<Integer,WorkspaceDataTable> openFile(String databaseName, long taxonomyID, int reactionID, double similarityThreshold) {

		List<String> alignCol = new ArrayList<>();
		alignCol.add("locus tag");
		alignCol.add("sequence id");
		alignCol.add("orthologous");
		alignCol.add("alignment score");
		alignCol.add("ec number");
		alignCol.add("query coverage");
		alignCol.add("target's coverage");
		alignCol.add("select");

		WorkspaceGenericDataTable data = new WorkspaceGenericDataTable(alignCol, "alignment table", "alignment results");

		Pair<Integer,WorkspaceDataTable> pair = new Pair<Integer, WorkspaceDataTable>(null, null);

		String path = FileUtils.getWorkspaceTaxonomyFillGapsFolderPath(databaseName, taxonomyID) + "FillGapReactions.txt";
		File file = new File(path);

		boolean pathExists = file.exists(); 
		if(pathExists){

			try {

				Scanner input = new Scanner(file);
				Map<String, String> credentials = LoadFromConf.loadReactionsThresholds(FileUtils.getConfFolderPath());
				similarityThreshold = Double.valueOf(credentials.get("similarity_threshold"));

				boolean empty = false;
				while (input.hasNextLine()){

					List<Object> line = new ArrayList<>();
					String[] str = input.nextLine().trim().split("\t");
					double simas = Double.parseDouble(str[1]);

					boolean existsReaction = Integer.parseInt(str[0])==reactionID;
					boolean hasSimilarity = simas<=similarityThreshold;
					boolean hasHomologue = !str[2].equals("empty");

					if(existsReaction && hasSimilarity ){
						if(hasHomologue) {
							for(int i=2; i < str.length; i++)
								line.add(str[i]);

							line.add(true);
							data.addLine(line);

						}
						else {
							empty = true;
						}
					}
				}

				pair.setB(data);

				input.close();
				if(data.getRowCount()!=0){
					pair.setA(1);
					return pair;

				}
				else if(data.getRowCount()==0 && empty==true) {
					pair.setA(0);
					return pair;
				}


			} catch (NumberFormatException e) {
				e.printStackTrace();
				pair.setA(-1);
				return pair;
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		pair.setA(-1);
		return pair;

	}



	public void initGui(){

		//JPanel jPanel = new JPanel();
		GridBagLayout thisLayout = new GridBagLayout();
		thisLayout.columnWeights = new double[] {0.0, 0.1, 0.0};
		thisLayout.columnWidths = new int[] {7, 7, 7};
		thisLayout.rowWeights = new double[] {0.0, 200.0, 0.0, 0.0, 0.0};
		thisLayout.rowHeights = new int[] {7, 50, 7, 3, 7};
		this.setLayout(thisLayout);
		this.setPreferredSize(new Dimension(900, 300));
		this.setSize(900, 300);
		this.setLayout(thisLayout);


		JScrollPane jScrollPane = new JScrollPane();
		this.add(jScrollPane,new GridBagConstraints(1,1, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

		JTable jTable = new JTable();
		jTable.setAutoCreateRowSorter(true);
		jTable.setModel(data);
		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment( JLabel.CENTER );
		jTable.setDefaultRenderer(String.class, centerRenderer);
		this.setTitle("select the row(s) to apply to the model");
		jScrollPane.setViewportView(jTable);
		this.setModal(true);


		JPanel jPanel = new JPanel();
		jPanel.setBorder(BorderFactory.createTitledBorder("options"));
		this.add(jPanel,new GridBagConstraints(1,3, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));



		jButtonOK = new JButton();
		jPanel.add(jButtonOK, new GridBagConstraints(2, 2, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		jButtonOK.setText("ok");
		jButtonOK.setToolTipText("apply to the model");
		jButtonOK.setIcon(new CreateImageIcon(new ImageIcon(getClass().getClassLoader().getResource("icons/Ok.png")),0.08).resizeImageIcon());
		jButtonOK.setPreferredSize(new Dimension(100, 30));
		jButtonOK.setSize(100, 30);
		jButtonOK.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {

				List<String> entries = readSelectedEntries(jTable);

				try {

					if(entries.size() > 0){
						apply(entries);
						simpleFinish();
					}
					else{
						windowToBack();
						Workbench.getInstance().warn("no rows selected to apply to the model!");
					}

				} catch (Exception e1) {
					windowToBack();
					Workbench.getInstance().warn("error while apllying!");
					e1.printStackTrace();
				}
			}
		});



		jButtonCancel = new JButton();
		jPanel.add(jButtonCancel, new GridBagConstraints(4, 2, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		jButtonCancel.setText("cancel");
		jButtonCancel.setToolTipText("close the tabel");
		jButtonCancel.setIcon(new CreateImageIcon(new ImageIcon(getClass().getClassLoader().getResource("icons/Cancel.png")),0.08).resizeImageIcon());
		jButtonCancel.setPreferredSize(new Dimension(100, 30));
		jButtonCancel.setSize(100, 30);
		jButtonCancel.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {

				simpleFinish();
			}
		});

	}

	/**
	 * Reads the selected lines.
	 * 
	 */
	private List<String> readSelectedEntries(JTable jTable){

		List<String> toApply =new ArrayList<>();

		for(int i = 0; i < jTable.getRowCount(); i++){

			Boolean selected = (boolean) jTable.getValueAt(i, 7);

			if( selected == true)
				toApply.add((String) jTable.getValueAt(i, 1));

		}
		return toApply;
	}

	private void windowToBack(){

		this.toBack();

	}

	private void simpleFinish() {
		this.setVisible(false);
		this.dispose();
	}
}