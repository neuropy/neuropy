# Code placed in load() in Experiments class to convert all ascii .din to binary of the same filename:

os.rename(self.path + self.name + '.din', self.path + self.name + '.din.txt')
fin = self.path + self.name + '.din.txt'
fout = self.path + self.name + '.din'
txtdin2binarydin(fin,fout)
os.remove(self.path + self.name + '.din.txt')
