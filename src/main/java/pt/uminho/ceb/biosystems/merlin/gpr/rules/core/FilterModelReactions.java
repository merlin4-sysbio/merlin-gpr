package pt.uminho.ceb.biosystems.merlin.gpr.rules.core;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.database.connector.databaseAPI.ModelAPI;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Connection;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.Enumerators.DatabaseType;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.H2DatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.database.connector.datatypes.MySQLDatabaseAccess;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.gpr.ReactionsGPR_CI;

/**
 * @author ODias
 *
 */
public class FilterModelReactions {

	private static final Logger logger = LoggerFactory.getLogger(FilterModelReactions.class);

	private DatabaseAccess dba;
	private DatabaseType database_type;
	private Map<String, Set<String>> databaseEnzymesReactions;
	private Map<String, String> annotations; 
	private boolean originalReactions;
	private Set<String> removed , kept, no_gpr;
	private Set<String> keptWithDifferentAnnotation;

	/**
	 * @param user
	 * @param password
	 * @param server
	 * @param port
	 * @param database
	 * @param originalReactions
	 */
	public FilterModelReactions(String user, String password, String server, int port, String database, boolean originalReactions) {
		
		if (this.database_type.equals(DatabaseType.MYSQL)) {
			this.dba = new MySQLDatabaseAccess(user, password, server, port, database);
		}else{
			this.dba = new H2DatabaseAccess(user, password, database, null);
		}

		new FilterModelReactions(dba, originalReactions);
	}

	/**
	 * @param dba
	 * @param originalReactions
	 */
	public FilterModelReactions(DatabaseAccess dba, boolean originalReactions) {

		this.dba = dba;
		this.databaseEnzymesReactions = new HashMap<>();
		this.originalReactions = originalReactions;

		this.removed = (new HashSet<>());
		this.kept = (new HashSet<>());
		this.keptWithDifferentAnnotation = (new HashSet<>());
		this.annotations  = new HashMap<>();
		this.no_gpr = (new HashSet<>());

	}


	/**
	 * @throws SQLException
	 */
	public Set<String> getReactionsFromModel() throws SQLException {

		Set<String> ret = new HashSet<String>();

		Connection conn = new Connection(this.dba);

		Statement stmt = conn.createStatement();

		ArrayList<String[]> result = ModelAPI.getReactionsFromModel(stmt, this.originalReactions);
		String [] list;
		
		for(int i=0; i<result.size(); i++) {	
			
			list = result.get(i);
			
			Set<String> reactions = new HashSet<>();

			if(this.databaseEnzymesReactions.containsKey(list[1]))
				reactions = this.databaseEnzymesReactions.get(list[1]);

			reactions.add(list[0]);
			ret.add(list[0]);

			this.databaseEnzymesReactions.put(list[1], reactions);
		}
		conn.closeConnection();
		return ret;
	}

	/**
	 * @param map
	 * @throws SQLException 
	 */
	public void filterReactions(Map<String, ReactionsGPR_CI> rpgs) throws SQLException {

		Map<String, Set<String>> gpr_map = new HashMap<>();
		Set<String> ecs = new HashSet<>();

		for(String reaction : rpgs.keySet()) {

			for(String key : rpgs.get(reaction).getProteins().keySet()) {

				ecs.add(key);

				if(rpgs.get(reaction).getProteins().get(key).isECnumberValid()) {

					this.annotations.put(reaction, rpgs.get(reaction).getProteins().get(key).getGeneRule());

					Set<String> reactions = new HashSet<>();

					if(gpr_map.containsKey(key))
						reactions = gpr_map.get(key);

					reactions.add(reaction);

					gpr_map.put(key, reactions);
				}
			}
		}

		Set<String> reactions = this.getReactionsFromModel();

		for(String ec : gpr_map.keySet()) {

			if(this.databaseEnzymesReactions.containsKey(ec)) {

				Set<String> r = this.databaseEnzymesReactions.get(ec);
				Set<String> r2 = new HashSet<>();

				for(String react : r) {

					for(String gpr_react : gpr_map.get(ec)) {

						if(react.contains(gpr_react)) {

							r2.add(react);
							kept.add(react);
							this.annotations.put(react, this.annotations.get(gpr_react));
						}
						else {

							removed.add(react);
						}
					}
				}
			}
			else {

				for(String gpr_react : gpr_map.get(ec)) {

					for(String react : reactions) {

						if(react.contains(gpr_react)) {

							keptWithDifferentAnnotation.add(react);
							this.annotations.put(react, this.annotations.get(gpr_react));
						}
					}
				}
			}
		}

		for(String ec : this.databaseEnzymesReactions.keySet()) {

			if(!gpr_map.containsKey(ec)) {

				if(this.databaseEnzymesReactions.containsKey(ec)) {

					Set<String> r = this.databaseEnzymesReactions.get(ec);

					for(String react : r) {

						this.no_gpr.add(react);
					}
				}
			}
		}

		this.keptWithDifferentAnnotation.removeAll(this.kept);
		this.removed.removeAll(this.kept);
		this.removed.removeAll(this.keptWithDifferentAnnotation);
		this.no_gpr.removeAll(this.kept);
		this.no_gpr.removeAll(this.keptWithDifferentAnnotation);
		this.removed.removeAll(this.no_gpr);

		String s = "Removed\t"+removed.size()+"\t"+removed+
				"\nKept\t"+kept.size()+"\t"+kept+
				"\nKept new annotation\t"+keptWithDifferentAnnotation.size()+"\t"+keptWithDifferentAnnotation+
				"\nNo GPR\t"+no_gpr.size()+"\t"+no_gpr;
		
		s+="\n\n";
		for(String r : this.annotations.keySet())
			s+="reaction:\t"+r+"\t"+this.annotations.get(r)+"\n";

		logger.info("Integration report: {}\nEnd Report.", s);

	}

	/**
	 * 
	 * @param keepReactionsWithNotes
	 * @param keepManualReactions
	 * @throws SQLException
	 */
	public void removeReactionsFromModel(boolean keepReactionsWithNotes, boolean keepManualReactions) throws SQLException {

		Map<String, String> notes_map = new HashMap<>();
		Set<String> reactionsToKeep = new HashSet<String>();
		Connection connection = new Connection(this.dba);
		Statement stmt = connection.createStatement();

		for (String name : this.removed) {

//			ResultSet rs = stmt.executeQuery("SELECT notes, isSpontaneous, isNonEnzymatic, source FROM reaction  WHERE reaction.name='"+name+"'");
			String[] result = ModelAPI.getReactionsInfo(stmt, name);
					
			if(result.length>0) {
				
				String old_note = result[0];

				if(old_note!=null && !old_note.isEmpty()) {

					if(keepReactionsWithNotes) {

						reactionsToKeep.add(name);
					}
					else {

						old_note = old_note.replaceAll(" \\| Removed by GPR rule","");
						old_note = old_note.replaceAll("Removed by GPR rule","");

						if(old_note.contains("new Annotation. automatic GPR")) {

							String[] data = old_note.split(" \\| ");

							if(data.length>2)
								old_note+=" | "+data[2];
						}

						if(old_note.contains("automatic GPR")) {

							String[] data = old_note.split(" \\| ");

							if(data.length>2)
								old_note+=" | "+data[2];
						}
						notes_map.put(name, old_note);
					}
				}
				if(Boolean.valueOf(result[1]) || Boolean.valueOf(result[2]))
					reactionsToKeep.add(name);

				if(keepManualReactions && result[3].equalsIgnoreCase("MANUAL"))
					reactionsToKeep.add(name);
			}
		}
		connection.closeConnection();

		logger.debug("Removed notes\t"+notes_map);

		java.sql.Connection conn = this.dba.openConnection();
		
		PreparedStatement statement = conn.prepareStatement("UPDATE reaction SET inModel=?, notes=? WHERE reaction.name=?");
		
		ModelAPI.removeReactionsFromModel(statement, this.removed, reactionsToKeep, notes_map);

		conn.close();
	}

	/**
	 * @throws SQLException
	 */
	public void setModelGPRsFromTool() throws SQLException {

		Connection connection = new Connection(this.dba);
		Statement stmt = connection.createStatement();

		Map<String, String> notes_map = ModelAPI.createNotesMap(stmt, this.kept);

		PreparedStatement statement = connection.prepareStatement("UPDATE reaction SET boolean_rule=?, notes=? WHERE reaction.name=?");

		ModelAPI.updateReactionTable(statement, this.kept, this.annotations, notes_map);

		ModelAPI.updateReactionTableWithDifferentAnnotation(statement, this.keptWithDifferentAnnotation, this.annotations, notes_map);
		
		connection.closeConnection();
	}
}
