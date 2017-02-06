class MySeq:

    def __init__(self,key,value):
        self.key = 'This is ' + key
        self.value = value.upper() 

    def __str__(self):
        rep = self.key + '\n' + self.value
        return rep
    #def get(self,key, value):
        #prin(self.key 

x = MySeq('human stuff','me! me! me!')
y = MySeq('feline stuff','catcatcatcat')  

print(x.key)
print(x.value,'\n')
print(y.key)
print(y.value,'\n')
print(x)
print(y)


