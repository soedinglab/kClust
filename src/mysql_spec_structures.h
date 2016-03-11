//sql_create_#(NAME, COMPCOUNT, SETCOUNT, TYPE1, ITEM1, ... TYPE#, ITEM#)
//compcount number of fields used for a comparison of two structures - not used
//setcount>0 automatically creates a constructor that takes setcount arguments and uses them for initialization of the first setcount elements in the structure 
sql_create_7(mysql_cluster, 1, 0,
    int, id,
    std::string, name,
	 int, size,
	 std::string, consensus_header,
    std::string, consensus_sequence,
    std::string, alignment,
    std::string, profil);

sql_create_4(mysql_sequence, 1, 0,
    int, id,
	 int, cluster,
    std::string, header,
    std::string, sequence);

