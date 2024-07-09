import psycopg2
from psycopg2 import sql
import random
from itertools import chain, combinations
import math
import numpy as np
from scipy import stats
from scipy.stats import skew

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
        print("Query executed successfully.")

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

def powerset(s):
    """
    Generates the powerset of a given set s.

    :param s: The input set
    :return: A list of subsets representing the powerset
    """
    s = list(s)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def generateGB(groupAtts):
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



def welch_ttest(sample1, sample2, alpha=0.05):
    """
    Perform Welch's t-test to determine if there is a significant difference
    between the means of two groups.

    Parameters:
    - sample1: list or numpy array of sample data for group 1
    - sample2: list or numpy array of sample data for group 2
    - alpha: significance level (default is 0.05)

    Returns:
    - t_stat: the calculated t-statistic
    - p_value: the two-tailed p-value
    - conclusion: string stating if we reject or fail to reject the null hypothesis
    """

    # Perform Welch's t-test
    t_stat, p_value = stats.ttest_ind(sample1, sample2, equal_var=False)

    # Conclusion
    if p_value < alpha:
        conclusion = "Reject the null hypothesis: There is a significant difference between the means of the two groups."
    else:
        conclusion = "Fail to reject the null hypothesis: There is no significant difference between the means of the two groups."

    return t_stat, p_value, conclusion



def compute_skewness(data):
    """
    Compute the skewness of a population.

    Parameters:
    - data: list or numpy array of sample data

    Returns:
    - skewness: the skewness of the data
    """

    # Convert data to a numpy array if it's not already
    #data = np.array(data)

    # Compute skewness
    skewness = skew(data, bias=False)

    return skewness

def claireStat(skew1, skew2, count1, count2, threshold=0.049):
    """
    Compute Claire statistics for testing if Welch test can be used

    Parameters:
    - skew1: skewness of first samples
    - skew2: skewness of second sample
    - count1: size of first sample
    - count2: size of second sample
    - threshold: why would we want to change?

    Returns:
    - skewness: boolean indicating if Welch test can be used
    """

    stat = abs(skew1/count1 - skew2/count2)
    if stat<threshold:
        return True
    else:
        return False




def ExampleUsage():

    # Sample data
    data = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5])
    # Compute skewness
    skewness = compute_skewness(data)
    # Print result
    print(f"Skewness of the data: {skewness}")

# Example usage

    # Sample data
    group1 = np.array([23, 45, 56, 67, 78, 89, 90, 56, 54, 67])
    group2 = np.array([35, 56, 76, 78, 45, 34, 23, 45, 67, 88])

    # Perform Welch's t-test
    t_stat, p_value, conclusion = welch_ttest(group1, group2)

    # Print results
    print(f"t-statistic: {t_stat}")
    print(f"p-value: {p_value}")
    print(conclusion)




#queryPattern="SELECT " + gb + "," + meas + " FROM " + table + " WHERE " + sel + " in " + vals + "group by " + gb +";"
#hypothesis=[('AA', 1),('UA', 2),('US', 3)]


table="fact_table"
measures=["nb_flights","departure_delay","late_aircraft"]
groupAtts=["departure_airport","date","departure_hour","flight","airline"]

measure="nb_flights"
selectionAtt="airline"
groupbyAtt=["departure_airport","date","departure_hour","flight"]
gb=""
sel="airline"
vals="('AA','UA','US')"
#meas="sum(nb_flights)"
meas="avg(departure_delay)"

sizeOfVals=5 #number of members for the hypothesis, the more the more likely hypothesis is ok?
q0=""
epsilon=0.05
alpha=0.05
p=0
H=[]
nbTests=0
threshold=0.1 # 10% of tuples violating the order
n=math.log(2/alpha,10) / pow(2,epsilon*epsilon)
print("n>= " + str(n))
n=math.ceil(n)

if __name__ == "__main__":

    # Database connection parameters
    dbname = "flight_dw"
    user = ""
    password = ""
    host = "localhost"
    port = "5432"

    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    if conn:
        #drawing a query
        pwsert = powerset(groupbyAtt)
        #l= len(pwsert)
        #print(pwsert)

        # empty group by set removed from powerset
        # since it is used to generate the hypothesis
        # note that hypothesis could also be user given
        pwsert.remove(())

        #compute top sizeOfVals members - could look for best combinations?
        queryEmptyGb=("SELECT " +  sel + ","
                 + " rank () over (  order by " + meas + " desc ) as rank" +
                 " FROM " + table +  " group by " + sel + " limit " + str(sizeOfVals) + ";")


        #queryEmptyGb=("SELECT " +  sel + ","
        #         + " rank () over (  order by " + meas + " desc ) as rank" +
        #         " FROM " + table + " WHERE " + sel + " in " + vals + " group by " + sel + ";")
        resultEmptyGb = execute_query(conn, queryEmptyGb)

        if resultEmptyGb is not None:
            print("Overall the ranking is as follows:")
            for row in resultEmptyGb:
                print(row)

        hypothesis=resultEmptyGb;
        vals=tuple([x[0] for x in hypothesis])
        #print("members: " + str(vals))

        for i in range(n):

            nb=random.randint(0, len(pwsert)-1)
            gb=pwsert[nb]
            strgb=""
            for i in range(len(gb)):
                strgb=strgb + str(gb[i])
                if i != len(gb)-1:
                    strgb=strgb+","

            print("group by is: " + strgb)

            # for debugging
            #strgb = "departure_airport"

            hyp = ""
            for i in range(len(hypothesis)):
                hyp = hyp + str(hypothesis[i])
                if i != len(hypothesis) - 1:
                    hyp = hyp + ","
            queryHyp = (
                    "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")

#        query = "SELECT " + strgb + "," + sel +"," + meas + " FROM " + table + " WHERE " + sel + " in " + vals + " group by " + strgb + "," + sel + " order by " + strgb + "," + sel + ";"
 #       query = ("SELECT " + strgb + "," + sel +"," + meas + ", "
  #               + " rank () over ( partition by " + strgb +  " order by " + meas + " desc )" +
   #              " FROM " + table + " WHERE " + sel + " in " + vals + " group by " + strgb + "," + sel +";")
        # EXCEPT ALL select * from (select strgb from fact_table) a  cross join (values hypotheses) as t ;
        #select * from (select departure_airport from fact_table) a  cross join (values ('AA',1),('AU',2)) as t ;

            query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
                 + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
                 " FROM " + table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ")


            queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "
                     +  table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ) x;")
            resultValues = execute_query(conn, queryValues)
            data=[]
            for row in resultValues:
                data.append(float(row[0]))
            nvalues=len(data)
            data=np.array(data)

            # Compute skewness
            skewness = compute_skewness(data)
            # Print result
            print(f"Skewness of the data: {skewness}")




            queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query +  " ) t3 except all " + queryHyp + ";")

            print("query is: " + queryExcept)

            queryCountGb=("select count(*) from (" + queryHyp + ") t4;")
            queryCountExcept=("select count(*) from (" + query + ") t5;")


        #strategy: use the db engine to check whether the hypothesis holds
        #cross join hypothesis to chosen group by set
        # compute the actual ranks of vals in the hypothesis for each group
        # except all the actual ranks with the hypothesis
        # remaining only the groups where hypothesis does not hold

        # Execute the query
            resultCountGb = execute_query(conn, queryCountGb)
            resultCountExcept = execute_query(conn, queryCountExcept)

            print("number of tuples checked: " + str( resultCountGb[0][0] ))
            print("number of exceptions: " + str(resultCountExcept[0][0]))
            print("ratio is: " + str(resultCountExcept[0][0]  / resultCountGb[0][0]))

            ratio=resultCountExcept[0][0]  / resultCountGb[0][0]

            # 1 if less than 10% errors
            # could keep actual ratio instead
            if ratio > threshold:
                H.append(0)
            else:
                H.append(1)

            nbTests=nbTests+1

            print("Expected value is: " + str(sum(H)/len(H)))
            print("Size of confidence interval around p: " + str(epsilon))
            print("Probability is of making a mistake: " + str(alpha))

            # Close the connection
        close_connection(conn)


