import psycopg2
from psycopg2 import sql

def connect_to_db(dbname, user, password, host='localhost', port='5432'):
    """
    Establishes a connection to the PostgreSQL database.

    :param dbname: Name of the database
    :param user: Database user
    :param password: User's password
    :param host: Database host address (default is 'localhost')
    :param port: Connection port number (default is '5432')
    :return: Connection object
    """
    try:
        conn = psycopg2.connect(
            dbname=dbname,
            user=user,
            password=password,
            host=host,
            port=port
        )
        print("Connection to database established successfully.")
        return conn
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return None


def execute_query(conn, query):
    """
    Executes a given SQL query using the established connection.

    :param conn: Connection object
    :param query: SQL query to be executed
    :return: Query result
    """
    try:
        cursor = conn.cursor()
        cursor.execute(query)
        conn.commit()
        #print("Query executed successfully.")

        try:
            result = cursor.fetchall()
            return result
        except psycopg2.ProgrammingError:
            # If the query does not return any data (like an INSERT or UPDATE)
            return None
    except Exception as e:
        print(f"Error executing query: {e}")
        return None
    finally:
        cursor.close()

def printResultSet(result):
    if result is not None:
        for row in result:
            print(row)

def close_connection(conn):
    """
    Closes the database connection.

    :param conn: Connection object
    """
    try:
        conn.close()
        print("Database connection closed.")
    except Exception as e:
        print(f"Error closing connection: {e}")



def generateGB(groupAtts, measures, table):
    """
    Generates group bys and sel from the list of all categorical attributes - unused so far

    :param groupAtts: all the categorical attributes
    :return:
    """
    for g in groupAtts:
        gb = groupAtts.remove(g)
        groupbyAtt = gb
        sel = g
        for m in measures:
            meas = "sum(" + m + ")"

            queryVals = ("select distinct " + sel + " from " + table + ";")


def getSample(conn, measBase, table, sel, sampleSize, method="SYSTEM", repeatable=False):
    # sampling using postgresql: https://www.postgresql.org/docs/current/sql-select.html#SQL-FROM
    # system (faster) is block based, bernouili (slower) is row based

    if repeatable:
        is_repeatable = 'REPEATABLE(42)'
    else:
        is_repeatable = ''

    querySample = (
            "SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + method + " (" + str(
        sampleSize) + ")" + is_repeatable + ";")
    # print(querySample)
    resultVals = execute_query(conn, querySample)
    print("Sample size in tuples: ", len(resultVals))
    return resultVals


def emptyGB(conn, nb_of_adom_vals, table, sel, meas):
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #               + " rank () over (  order by " + meas + " desc ) as rank" +
    #              " FROM " + table +  " WHERE " + sel + " in (" + hyp + ") group by " + sel +  ";")

    #hyp = ""
    #for i in range(len(vals)):
    #    hyp = hyp + "'" + str(vals[i]) + "'"
    #    if i != len(vals) - 1:
    #        hyp = hyp + ","

    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel + " limit " + str(nb_of_adom_vals) + ";")

    #print(queryEmptyGb)
    resultEmptyGb = execute_query(conn, queryEmptyGb)

    queryEmptyGbAll = ("SELECT " + sel + ","
                       + " rank () over (  order by " + meas + " desc ) as rank" +
                       " FROM " + table + " group by " + sel + ";")
    resultEmptyGbAll = execute_query(conn, queryEmptyGbAll)

    return resultEmptyGb, resultEmptyGbAll