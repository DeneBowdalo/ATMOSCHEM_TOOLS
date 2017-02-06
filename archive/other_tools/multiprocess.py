from multiprocessing import Process

def say_hello(name='world'):
    print "Hello, %s" % name
names = ['john', 'dave', 'bob']

p = Process(target=say_hello,args= names)
p.start()
p.join()
