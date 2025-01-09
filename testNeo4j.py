import configparser
import json
from neo4j import GraphDatabase


class Neo4jConnection:
    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()

    def query(self, query, parameters=None):
        with self._driver.session() as session:
            result = session.run(query, parameters)
            return [record.data() for record in result]



if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read('configs/flightsNeo4j.ini')


    USER = "PM"

    # Database connection parameters
    dbname = config[USER]['dbname']
    user = config[USER]['user']
    password = config[USER]['password']
    host = config[USER]['host']
    port = int(config[USER]['port'])

    # Cube info
    table = config["Common"]['table']
    measures = json.loads(config.get("Common", "measures"))
    groupbyAtt = json.loads(config.get("Common", "groupbyAtt"))
    sel = config["Common"]['sel']
    meas = config["Common"]['meas']
    measBase = config["Common"]['measBase']
    function = config["Common"]['function']
    prefs = json.loads(config.get("Common", "preferred"))

    # Replace these with your Neo4j instance details
    #NEO4J_URI = "bolt://localhost:7687"  # Typically starts with "bolt://"
    NEO4J_URI = "bolt://"+str(host)+":"+ str(port)  # Typically starts with "bolt://"


    # Example query
    QUERY = """
    MATCH (n)
    RETURN n LIMIT 5
    """

    # Initialize the connection
    conn = Neo4jConnection(uri=NEO4J_URI, user=user, password=password)

    try:
        # Execute the query and print results
        results = conn.query(QUERY)
        print("Query Results:")
        for result in results:
            print(result)
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Close the connection
        conn.close()


