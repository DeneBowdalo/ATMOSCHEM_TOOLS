class NotaList(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Vector(object):
    def setdata(self,vector_name=[]):
        self.vector_name = vector_name
        if type(vector_name) != list:
            raise NotaList('Input is not a list')
        for x in self.vector_name:
           if type(x) != int:
               raise NotaList('List needs to be numbers')
   
    def getdata(self):
        print(self.vector_name)

    def add(self,second_vector):
        added_vectors = []
        for i in range(len(self.vector_name)):
            added_vectors.append(self.vector_name[i] + second_vector.vector_name[i])
        added_object = Vector()
        added_object.setdata(added_vectors)
        return added_object


vec1 = Vector()
vec2 = Vector()

vec1.setdata([1,2,3])
vec2.setdata([2,3,4])

vec1.getdata()
vec2.getdata()

vec3 = vec1.add(vec2)
vec3.getdata()
