wList = ['c','d']
word = 'a_b_c'
word.split("_")
print word.split("_")[0]
#word = word[len(word.split("_")[0])+1:]
print word
if not any( w in word for w in wList): print 'test'
