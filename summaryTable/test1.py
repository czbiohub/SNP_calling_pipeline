"""defining the sclaud class!!"""
class sclaud_obj:

	type = 'sclaud_obj'

    # defining a function inherited by sclaud_obj
	def f(self):
		return "you've just created an object of the class sclaud_obj"

	def __init__(self, name):
		self.data = []
		self.name = name

x = sclaud_obj('myLAUD_obj')

print(x.f())
print(x.type)
print(x.name)