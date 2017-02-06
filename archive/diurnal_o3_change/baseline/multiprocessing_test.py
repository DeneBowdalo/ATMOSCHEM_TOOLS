import multiprocessing
from multiprocessing import Queue
import matplotlib.pyplot as plt
from functools import partial
import numpy as np

def g(tup):
    return simple(*tup)

def simple(n,a,b):
    #print 'a =', a
    #print 'b =', b
    #print a
    #print b
    return (a,b,n)


a = np.arange(0,10000000,1)
a = np.vstack((a,a,a,a,a))

#print a.shape

b = np.arange(10000000,20000000,1)
b = np.vstack((b,b,b,b,b))
#b = np.vstack((b,b))

print b.shape

i = np.vstack((a,b))

#print i.shape
#print i[0]

#i = i.tolist()


print i.shape

#partial_simple = partial(simple, a, b)

if __name__ == '__main__':
    #queue1 = Queue()
    #for x in range(5):
    pool = multiprocessing.Pool(processes=6)
    #pool_outputs = pool.map(simple, inputs)
    
    #p = multiprocessing.Process(target=simple, args=(a[x],b[x],queue1))
    #print p
    #p.start()
        #ss =  queue1.get()
        #print multiprocessing.active_children()
        #try:
        #    big_array = np.vstack((big_array,ss))
        #except:
        #    big_array = ss
        #p.join()
    
    #big_array = np.array(big_array)
    #print big_array.shape
    results = [pool.apply_async(simple, (x,a[x],b[x])) for x in range(5)]
    roots = [r.get() for r in results]
    print roots
    #print roots
    #ans = pool.map(lambda args: simple(*args), i)
    #ans = pool.map(partial(simple, b), a)
    #pool.close()
    #pool.join()
    roots = np.array(roots)
    print roots.shape
    print roots[1,2]
    #print roots[0]
