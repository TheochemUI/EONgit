#!/usr/bin/env python

f = open('barriers.test')
ref_barriers = []
for line in f:
    ref_barriers.append(float(line))
f.close()

f = open('states/0/processtable')
my_barriers = []
for line in f:
    fields = line.split()
    if fields[0] == 'proc': continue
    barrier = round(float(fields[6]),3)
    my_barriers.append(barrier)

for ref_barrier in ref_barriers:
    if ref_barrier not in my_barriers:
        print "missing barrier %.3f" % ref_barrier
