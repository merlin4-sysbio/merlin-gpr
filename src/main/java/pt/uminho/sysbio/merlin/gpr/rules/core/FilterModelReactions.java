package pt.uminho.sysbio.merlin.gpr.rules.core;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.sysbio.common.database.connector.databaseAPI.containers.gpr.ReactionsGPR_CI;
import pt.uminho.sysbio.common.database.connector.datatypes.Connection;
import pt.uminho.sysbio.common.database.connector.datatypes.DatabaseAccess;
import pt.uminho.sysbio.common.database.connector.datatypes.H2DatabaseAccess;
import pt.uminho.sysbio.common.database.connector.datatypes.MySQLDatabaseAccess;
import pt.uminho.sysbio.common.database.connector.datatypes.Enumerators.DatabaseType;

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
			this.dba = new H2DatabaseAccess(user, password, database);
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

		ResultSet rs = stmt.executeQuery("SELECT * FROM reaction " +
				"INNER JOIN reaction_has_enzyme  ON (reaction_idreaction = reaction.idreaction) " +
				" WHERE reaction.inModel AND originalReaction="+this.originalReactions
				);

		while (rs.next()) {

			Set<String> reactions = new HashSet<>();

			if(this.databaseEnzymesReactions.containsKey(rs.getString("enzyme_ecnumber")))
				reactions = this.databaseEnzymesReactions.get(rs.getString("enzyme_ecnumber"));

			reactions.add(rs.getString("name"));
			ret.add(rs.getString("name"));

			this.databaseEnzymesReactions.put(rs.getString("enzyme_ecnumber"), reactions);
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
			s+=r+"\t"+this.annotations.get(r)+"\n";

		logger.info("{}", s);

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

			ResultSet rs = stmt.executeQuery("SELECT notes, isSpontaneous, isNonEnzymatic, source FROM reaction  WHERE reaction.name='"+name+"'");

			if(rs.next()) {

				String old_note = rs.getString(1);

				if(old_note!=null && !old_note.isEmpty()) {

					if(keepReactionsWithNotes) {

						reactionsToKeep.add(name);
					}
					else {

						old_note = old_note.replaceAll(" \\| Removed by GPR tool","");
						old_note = old_note.replaceAll("Removed by GPR tool","");

						if(old_note.contains("New Annotation. GPR set from tool")) {

							String[] data = old_note.split(" \\| ");

							if(data.length>2)
								old_note+=" | "+data[2];
						}

						if(old_note.contains("GPR set from tool")) {

							String[] data = old_note.split(" \\| ");

							if(data.length>2)
								old_note+=" | "+data[2];
						}
						notes_map.put(name, old_note);
					}
				}


				if(rs.getBoolean(2) || rs.getBoolean(3))
					reactionsToKeep.add(name);

				if(keepManualReactions && rs.getString(4).equalsIgnoreCase("MANUAL"))
					reactionsToKeep.add(name);
			}
		}
		connection.closeConnection();

		logger.debug("Removed notes\t"+notes_map);

		java.sql.Connection conn = this.dba.openConnection();
		PreparedStatement statement = conn.prepareStatement("UPDATE reaction SET inModel=?, notes=? WHERE reaction.name=?");

		int i = 0;
		for (String name : this.removed) {

			if(!reactionsToKeep.contains(name)) {

				String note = "";

				if(notes_map.containsKey(name)) {

					note = notes_map.get(name)+ " | ";
				}

				note += "Removed by GPR tool";

				statement.setString(1, "false");
				statement.setString(2, note);
				statement.setString(3, name);

				statement.addBatch();

				if ((i + 1) % 1000 == 0) {

					statement.executeBatch(); // Execute every 1000 items.
				}
				i++;
			}
		}
		statement.executeBatch();

		conn.close();
	}

	/**
	 * @throws SQLException
	 */
	public void setModelGPRsFromTool() throws SQLException {

		Map<String, String> notes_map = new HashMap<>();
		Map<String, String> notes_map_new = new HashMap<>();

		Connection connection = new Connection(this.dba);
		Statement stmt = connection.createStatement();

		for (String name : this.kept) {

			ResultSet rs = stmt.executeQuery("SELECT notes FROM reaction WHERE reaction.name='"+name+"'");

			if(rs.next() && rs.getString(1)!=null && !rs.getString(1).isEmpty())
				notes_map.put(name, rs.getString(1));
		}

		for (String name : this.keptWithDifferentAnnotation) {

			ResultSet rs = stmt.executeQuery("SELECT notes FROM reaction  WHERE reaction.name='"+name+"'");

			if(rs.next() && rs.getString(1)!=null && !rs.getString(1).isEmpty())
				notes_map_new.put(name, rs.getString(1));
		}

		connection.closeConnection();

		java.sql.Connection conn = this.dba.openConnection();

		PreparedStatement statement = conn.prepareStatement("UPDATE reaction SET boolean_rule=?, notes=? WHERE reaction.name=?");

		int i = 0;
		for (String name : this.kept) {

			//String note = "GENE_ASSOCIATION: " + this.annotations.get(name)+" | GPR set from tool";
			String note = this.annotations.get(name);

			String old_note = "GPR set from tool";
			
			if(notes_map.containsKey(name)) {
				
				old_note = notes_map.get(name);

				if(!old_note.contains("GPR set from tool"))
					old_note = old_note.trim().concat(" | GPR set from tool");
			}

			statement.setString(1, note);
			statement.setString(2, old_note);
			statement.setString(3, name);
			statement.addBatch();

			if ((i + 1) % 1000 == 0) {

				statement.executeBatch(); // Execute every 1000 items.
			}
			i++;
		}
		statement.executeBatch();

		i = 0;
		for (String name : this.keptWithDifferentAnnotation) {

			//String note = "GENE_ASSOCIATION: " + this.annotations.get(name)+" | New Annotation. GPR set from tool.";
			String note = this.annotations.get(name);

			String old_note = "New Annotation. GPR set from tool";
			
			if(notes_map.containsKey(name)) {
				
				old_note = notes_map.get(name);

				if(!old_note.contains("New Annotation. GPR set from tool"))
					old_note = old_note.trim().concat(" | GPR set from tool");
			}

			statement.setString(1, note);
			statement.setString(2, old_note);
			statement.setString(3, name);
			statement.addBatch();

			if ((i + 1) % 1000 == 0) {

				statement.executeBatch(); // Execute every 1000 items.
			}
			i++;
		}
		statement.executeBatch();

		conn.close();
	}
}
