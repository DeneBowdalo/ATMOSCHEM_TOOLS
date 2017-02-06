class Vector:
    def setdata(self, **vector_list):
        self.data = vector_list
    def getdata(self):
        return self.data 
    def add(self, **other_self):
        other_self = self.data + self.data

