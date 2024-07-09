import psycopg2
from psycopg2 import sql
import random
from itertools import chain, combinations
import math
from scipy import stats
from scipy.stats import skew
import numpy as np
from scipy.stats import f

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
        # test can be used
        return True
    else:
        return False





def brown_forsythe_test(*groups):
    """
    Perform the Brown-Forsythe test for equal variances.

    Parameters:
    *groups : array_like
        Arrays of sample data. There must be at least two arrays.

    Returns:
    W : float
        The Brown-Forsythe test statistic.
    p_value : float
        The p-value for the test.
    """
    # Number of groups
    k = len(groups)

    if k < 2:
        raise ValueError("At least two groups are required")

    # Calculate the group sizes
    n = np.array([len(group) for group in groups])
    N = np.sum(n)

    # Calculate the group medians
    medians = np.array([np.median(group) for group in groups])

    # Calculate the absolute deviations from the medians
    Z = [np.abs(group - median) for group, median in zip(groups, medians)]

    # Calculate the overall mean of the absolute deviations
    Z_flat = np.concatenate(Z)
    Z_mean = np.mean(Z_flat)

    # Calculate the between-group sum of squares
    SSB = np.sum(n * (np.array([np.mean(z) for z in Z]) - Z_mean) ** 2)

    # Calculate the within-group sum of squares
    SSW = np.sum([np.sum((z - np.mean(z)) ** 2) for z in Z])

    # Calculate the Brown-Forsythe test statistic
    dfb = k - 1
    dfw = N - k
    W = (SSB / dfb) / (SSW / dfw)

    # Calculate the p-value
    p_value = f.sf(W, dfb, dfw)

    return W, p_value





def ExampleUsage():

# Example usage for skewness

    # Sample data
    data = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5])
    # Compute skewness
    skewness = compute_skewness(data)
    # Print result
    print(f"Skewness of the data: {skewness}")

# Example usage for Welch

    # Sample data
    group1 = np.array([23, 45, 56, 67, 78, 89, 90, 56, 54, 67])
    group2 = np.array([35, 56, 76, 78, 45, 34, 23, 45, 67, 88])

    # Perform Welch's t-test
    t_stat, p_value, conclusion = welch_ttest(group1, group2)

    # Print results
    print(f"t-statistic: {t_stat}")
    print(f"p-value: {p_value}")
    print(conclusion)

# Example usage for Brown Forsythe
    group1 = [23, 20, 25, 22, 21, 20, 23, 24]
    group2 = [22, 21, 24, 21, 20, 19, 21, 23]
    group3 = [26, 24, 27, 25, 26, 27, 28, 29]

    W, p_value = brown_forsythe_test(group1, group2, group3)
    print(f"Brown-Forsythe test statistic: {W}")
    print(f"P-value: {p_value}")


def generateHypothesis(conn):
    # compute top sizeOfVals members - could look for best combinations?
    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")

    resultEmptyGb = execute_query(conn, queryEmptyGb)

    """
    if resultEmptyGb is not None:
        print("Overall the ranking is as follows:")
        for row in resultEmptyGb:
            print(row)
    """

    return resultEmptyGb

def generateRandomQuery(pwsert,hypothesis):
    nb = random.randint(0, len(pwsert) - 1)
    gb = pwsert[nb]
    strgb = ""
    for i in range(len(gb)):
        strgb = strgb + str(gb[i])
        if i != len(gb) - 1:
            strgb = strgb + ","

    print("group by is: " + strgb)

    # for debugging
    # strgb = "departure_airport"

    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    queryHyp = (
            "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ")

    queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "
                   + table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ) x;")


    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")

    #return query, queryHyp, queryValues, queryExcept, strgb
    return queryValues,queryCountGb,queryCountExcept


def getValues(queryValues, vals, v, conn):
    queryVal = queryValues.replace(str(vals), "('" + v + "')")
    resultValues = execute_query(conn, queryVal)
    data = []
    for row in resultValues:
        data.append(float(row[0]))

    np.array(data)
    return data


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

sizeOfVals=5 #number of members for the hypothesis
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
        #compute powerset of categorical attributes
        pwsert = powerset(groupbyAtt)

        # generating the hypothesis
        # hypothesis is a ranking of members

        # empty group by set removed from powerset
        # since it is used to generate the hypothesis
        # TODO hypothesis could also be user given
        pwsert.remove(())

        hypothesis=generateHypothesis(conn);
        vals=tuple([x[0] for x in hypothesis])

        # n queries enough according to Hoeffding
        for i in range(n):

            # generate the random query
            #query, queryHyp, queryValues, queryExcept, strgb = generateRandomQuery(pwsert,hypothesis)
            queryValues,queryCountGb,queryCountExcept=generateRandomQuery(pwsert,hypothesis)


            # some statistics on the random query result: skew and size
            resultValues = execute_query(conn, queryValues)
            data=[]
            for row in resultValues:
                data.append(float(row[0]))
            nvalues=len(data)
            data=np.array(data)
            # Compute skewness
            skewness = compute_skewness(data)

            print("Size of the data: " + str(nvalues))
            print(f"Skewness of the data: {skewness}")

            S=[]
            # compute stats for every member of the hypothesis
            for v in vals:
                queryVal=queryValues.replace(str(vals),"('" + str(v) + "')")
                #print(queryVal)
                resultValues = execute_query(conn, queryVal)
                data = []
                for row in resultValues:
                    data.append(float(row[0]))
                nvalues = len(data)
                data = np.array(data)
                # Compute skewness
                skewness = compute_skewness(data)
                # Print result
                #print("Size of the data: " + str(nvalues))
                #print(f"Skewness of the data: {skewness}")
                S.append((v,nvalues,skewness))

            #print(S)
            # example of statstical test for first 2 members
            b=claireStat(S[0][2], S[1][2], S[0][1], S[1][1])
            print("for " + S[0][0] + " and " + S[1][0] + " Claire test says: " + str(b))
            if b:
                print("Welch test can be used")
                sample1 = getValues(queryValues,vals, S[0][0], conn)
                sample2 = getValues(queryValues,vals, S[1][0], conn)
                t_stat, p_value, conclusion=welch_ttest(sample1, sample2)
                print(conclusion)


            #strategy: use the db engine to check whether the hypothesis holds
            #cross join hypothesis to chosen group by set
            # compute the actual ranks of vals in the hypothesis for each group
            # except all the actual ranks with the hypothesis
            # remaining only the groups where hypothesis does not hold


            resultCountGb = execute_query(conn, queryCountGb)
            resultCountExcept = execute_query(conn, queryCountExcept)

            print("number of tuples checked: " + str( resultCountGb[0][0] ))
            print("number of exceptions: " + str(resultCountExcept[0][0]))
            print("ratio is: " + str(resultCountExcept[0][0]  / resultCountGb[0][0]))

            ratio=resultCountExcept[0][0]  / resultCountGb[0][0]


            # keep actual ratio
            # could use Bernouilli random variables instead: 1 if less than 10% errors
            if ratio > threshold:
                H.append(ratio)
            #    H.append(0)
            else:
                H.append(ratio)
            #    H.append(1)

            nbTests=nbTests+1

            print("Expected value is: " + str(sum(H)/len(H)))
            print("Size of confidence interval around p: " + str(epsilon))
            print("Probability is of making a mistake: " + str(alpha))

            # Close the connection
        close_connection(conn)


