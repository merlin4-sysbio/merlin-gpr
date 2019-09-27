package pt.uminho.ceb.biosystems.merlin.gpr.rules.core;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.core.containers.gpr.ReactionsGPR_CI;
import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.ModelAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Connection;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelReactionsServices;

/**
 * @author ODias
 *
 */
public class FilterModelReactions {

	private static final Logger logger = LoggerFactory.getLogger(FilterModelReactions.class);

	private Connection connection;
	private Map<String, Set<String>> databaseEnzymesReactions;
	private Map<String, String> annotations; 
	private boolean originalReactions;
	private Set<String> removed , kept, no_gpr;
	private Set<String> keptWithDifferentAnnotation;

	private String workspaceName;

	/**
	 * @param dba
	 * @param originalReactions
	 */
	public FilterModelReactions(Connection connection, String workspaceName, boolean originalReactions) {

		this.connection = connection;
		this.workspaceName = workspaceName;
		this.databaseEnzymesReactions = new HashMap<>();
		this.originalReactions = originalReactions;

		this.removed = (new HashSet<>());
		this.kept = (new HashSet<>());
		this.keptWithDifferentAnnotation = (new HashSet<>());
		this.annotations  = new HashMap<>();
		this.no_gpr = (new HashSet<>());

	}

	/**
	 * 
	 * @param keepReactionsWithNotes
	 * @param keepManualReactions
	 * @throws SQLException
	 */
	public void removeReactionsFromModel(boolean keepReactionsWithNotes, boolean keepManualReactions) throws Exception {

		ModelReactionsServices.removeReactionsFromModelByBooleanRule(this.workspaceName, 
				new ArrayList<String>(this.removed), keepReactionsWithNotes, keepManualReactions);
	}

	/**
	 * @throws SQLException
	 */
	public void setModelGPRsFromTool() throws Exception {
		
		ModelReactionsServices.updateBooleanRuleAndNotes(this.workspaceName, new ArrayList<String>(this.kept),
				this.annotations, "automatic GPR");
		
		ModelReactionsServices.updateBooleanRuleAndNotes(this.workspaceName, new ArrayList<String>(this.keptWithDifferentAnnotation), 
				this.annotations, "New Annotation. GPR set from tool");

	}
}
