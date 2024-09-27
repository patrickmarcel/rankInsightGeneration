import math

# delta is the probabiliy of making an error, sigma observed variance and t error
def sizeOfSample(delta, sigma, t):
    return (math.log(2/delta) * (2* math.pow(sigma,2) + (2*t)/3)) / math.pow(t,2)


# probability of making a mistake for error t
def bernsteinBound(delta, sigma, t):
    return 2* math.exp(-(math.pow(t,2) /2)/(math.pow(sigma,2) +t/3))

#delta the probability of making an error
def bersteinError(delta, sigma):
    return math.sqrt(math.ln(2/delta)*sigma/2) + (math.ln(2/delta)*(1+math.sqrt(2)))/math.sqrt(2)*6