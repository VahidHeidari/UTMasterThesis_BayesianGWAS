import time



IS_VERBOSE         = True



def Log(msg, is_verbose=IS_VERBOSE):
	if not IS_VERBOSE:
		return

	time_str = time.ctime(time.time())
	out_msg = time_str + ' ' + str(msg)
	print(out_msg)
	with open('log.txt', 'a') as f:
		f.write(out_msg + '\n')
		f.flush()


def LogError(msg):
	Log(msg, True)



if __name__ == '__main__':
	print('This is a module!')
