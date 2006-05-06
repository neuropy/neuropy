# This code was placed in load() in Experiments class to rename all experiment files to 0-based ids
# re came in handy!

basename = self.name # doesn't include any .din or .textheader or .sweeptable, etc. extensions
filelist = os.listdir(self.path)
p = re.compile('\.[0-9]+\.') # this pattern is any number of digits between a pair of periods
for fname in filelist:
	if fname.startswith(basename):
		try:
			(si,ei) = p.search(fname).span() # search for the pattern and return the start and end indices of the match in a tuple
			si += 1 # one past the first period in the match
			ei -= 1 # one before the last period in the match
			expid = int(fname[si:ei]) # gets you the 1-based expid
			expid -= 1 # convert to 0-based
			ndigits = ei - si
			expid = str(expid).zfill(ndigits) # convert back to string and pad with leading zeros if appropriate
			newfname = fname[0:si] + expid + fname[ei:] # this includes any extensions on the filename
			print fname, '  >>  ', newfname
			os.rename(self.path+fname, self.path+newfname)
		except AttributeError: # there was no match on this fname, search returned None, can't perform .span() on None
			pass
