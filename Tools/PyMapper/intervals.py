
def HasIntersection(interval_1, interval_2):
	return interval_1[0] <= interval_2[1] and interval_2[0] <= interval_1[1]


def IsInside(interval_small, interval_big):
	return interval_big[0] <= interval_small[0] and interval_small[1] <= interval_big[1]



class Intervals:
	def __init__(self):
		self.intervals = []


	def AddInterval(self, st, ed):
		if ed < st:
			tmp = st; st = ed; ed = st

		intervals_len = len(self.intervals)
		if intervals_len == 0:
			self.intervals.append((st, ed))
			return

		new_interval = (st, ed)
		for i in range(intervals_len):
			if HasIntersection(self.intervals[i], new_interval):
				if IsInside(self.intervals[i], new_interval):
					return

				self.Merge(i, new_interval)
				while i + 1 < intervals_len and HasIntersection(self.intervals[i + 1], self.intervals[i]):
					self.Merge(i, self.intervals[i + 1])
					self.intervals.pop(i + 1)
				return
			elif new_interval[1] < self.intervals[i][0]:
				self.intervals.insert(i, new_interval)
				return

		self.intervals.append(new_interval)


	def Merge(self, idx, new_interval):
		st = self.intervals[idx][0]
		ed = self.intervals[idx][1]
		self.intervals[idx] = (min(st, new_interval[0]), max(ed, new_interval[1]))


if __name__ == '__main__':
	print('This is a module!')
