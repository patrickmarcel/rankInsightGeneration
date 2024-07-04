import random
from itertools import chain, combinations
import math

from db_utils import connect_to_db, execute_query, close_connection


def powerset(s):
    """
    Generates the powerset of a given set s.

    :param s: The input set
    :return: A list of subsets representing the powerset
    """
    s = list(s)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


table = "fact_table"
measure = "nb_flights"
selectionAtt = "airline"
groupbyAtt = ["departure_airport", "date", "departure_hour", "flight"]
gb = ""
sel = "airline"
vals = "('AA','UA','US')"
meas = "sum(nb_flights)"
#queryPattern="SELECT " + gb + "," + meas + " FROM " + table + " WHERE " + sel + " in " + vals + "group by " + gb +";"

#hypothesis=[('AA', 1),('UA', 2),('US', 3)]
q0 = ""
epsilon = 0.05
alpha = 0.05
p = 0
H = []
nbTests = 0
threshold = 0.1
n = math.log(2 / alpha, 10) / pow(2, epsilon * epsilon)
print("n>= " + str(n))
n = math.ceil(n)


def get_hypothesis():
    global hypothesis
    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " WHERE " + sel + " in " + vals + " group by " + sel + ";")
    resultEmptyGb = execute_query(conn, queryEmptyGb)
    if resultEmptyGb is not None:
        print("Overall the ranking is as follows:")
        for row in resultEmptyGb:
            print(row)
    return resultEmptyGb


if __name__ == "__main__":
    # Database connection parameters
    dbname = "covid"
    user = "marcel"
    password = "1598marcel"
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
        print("len of power set", len(pwsert))

        hypothesis = get_hypothesis()

        for i in range(n):

            nb = random.randint(0, len(pwsert) - 1)
            gb = pwsert[nb]
            strgb = ",".join(gb)

            # for debugging
            #strgb = "departure_airport"

            hyp = ",".join(map(str, hypothesis))
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
                     " FROM " + table + " WHERE " + sel + " in " + vals + " group by " + strgb + "," + sel + " ")

            queryExcept = (
                        "select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + ";")

            print("query is: " + queryExcept)

            queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
            queryCountExcept = ("select count(*) from (" + query + ") t5;")

            # strategy: use the query engine to check whether the hypothesis holds
            # cross join hypothesis to chosen group by set
            # compute the actual ranks of vals in the hypothesis for each group
            # except all the actual ranks with the hypothesis
            # remaining only the groups where hypothesis does not hold

            # Execute the query
            resultCountGb = execute_query(conn, queryCountGb)
            resultCountExcept = execute_query(conn, queryCountExcept)

            ratio = resultCountExcept[0][0] / resultCountGb[0][0]

            print("number of tuples checked: " + str(resultCountGb[0][0]))
            print("number of exceptions: " + str(resultCountExcept[0][0]))
            print("ratio is: " + str(ratio))

            if ratio > threshold:
                H.append(0)
            else:
                H.append(1)

            nbTests = nbTests + 1

            #       if result is not None:
            #           for row in result:
            #               print(row)

            print("probability is: " + str(sum(H) / len(H)))

            # Close the connection
        close_connection(conn)
