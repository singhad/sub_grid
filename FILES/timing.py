#python3
import atexit
from time import time, strftime, localtime
from datetime import timedelta

def secondsToStr(elapsed=None):
    if elapsed is None:
        return strftime("%Y-%m-%d %H:%M:%S", localtime())
    else:
        return str(timedelta(seconds=elapsed))

def log(s, elapsed=None):
    line = "="*40
    print(line)
    print(secondsToStr(), '-', s)
    if elapsed:
        print("Elapsed time:", elapsed)
    print(line)
    print()

def endlog():
    end = time()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))

start = time()
atexit.register(endlog)
log("Start Program")

# python2
# import atexit
# from time import clock
#
# def secondsToStr(t):
#     return "%d:%02d:%02d.%03d" % \
#         reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
#             [(t*1000,),1000,60,60])
#
# line = "="*40
# def log(s, elapsed=None):
#     print line
#     print secondsToStr(clock()), '-', s
#     if elapsed:
#         print "Elapsed time:", elapsed
#     print line
#     print
#
# def endlog():
#     end = clock()
#     elapsed = end-start
#     log("End Program", secondsToStr(elapsed))
#
# def now():
#     return secondsToStr(clock())
#
# start = clock()
# atexit.register(endlog)
# log("Start Program")
