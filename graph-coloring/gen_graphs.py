
import random
import sys

def generate_ring(nodes):
	#print nodes   #nodes
	#print nodes   #edges
	for i in range(nodes):
		if i != nodes - 1:
			print (i, i + 1)
		if i == 0:
			print (i, nodes - 1)
	
					
def generate_crown(halfnodes):
	print (2 * halfnodes)
	print (halfnodes * (halfnodes - 1))
	for i in range(halfnodes):
		a = 2 * i
		for j in range(halfnodes):
			b = 2 * j + 1
			if (i != j):
				print (a, b)
				
				
def generate_random_graph(nodes, p):
	list_edges = []
	for i in range(nodes):
		for j in range(i + 1, nodes, 1):
			rand_num = random.random()
			if rand_num < p:
				list_edges.append(str(i) + " " + str(j))
	print (nodes)
	print (len(list_edges))
	for s in list_edges:
		print (s)
			

# m - x
# n - y
# p - z
def generate_grid(m, n, p):
	print (n * m * p)
	edges = p * ((m - 1) * n + (n - 1) * m) + (p - 1) * m * n
	print (edges)
	num = 0
	for z in range(0, p, 1):
		for y in range(0, n, 1):
			for x in range(0, m, 1):
				id = z * (m * n) + y * m + x
				if x < m - 1:
					print (id, id + 1)
					num += 1
				if y < n - 1:
					print (id, id + m)
					num += 1
				if z < p - 1:
					print (id, id + (m * n))
					num += 1
	if num != edges:
		sys.stderr.write("error in edges")
		


def little_test():
	generate_ring(40)
	#print
	generate_crown(3)
	#print
	generate_random_graph(10, 0.3)
	#print
	generate_cube(30, 40, 20)


def input_err():
	sys.stderr.write("Parameter Error\n")
	sys.exit(1)

if len(sys.argv) < 2:
	input_err()

graph = sys.argv[1]

if graph == "ring":
	if len(sys.argv) != 3:
		input_err()
	tt = int(sys.argv[2])
	generate_ring(tt)
elif graph == "crown":	
	if len(sys.argv) != 3:
		input_err()
	meiott = int(sys.argv[2])
	generate_crown(meiott)
elif graph == "random":	
	if len(sys.argv) != 4:
		input_err()
	tt = int(sys.argv[2])
	probb = float(sys.argv[3])
	generate_random_graph(tt, probb)
elif graph == "grid":
	if len(sys.argv) != 5:
		input_err()
	ttx = int(sys.argv[2])
	tty = int(sys.argv[3])
	ttz = int(sys.argv[4])
	generate_grid(ttx, tty, ttz)
else:
	input_err()
	