#!/bin/python

species_list = open("bacterial.txt","r")
uniq_species = open("uniq.sp.txt","w")
species_kinds = []
with species_list as fh:
	for line in fh:
		data = line.strip().split("_")
		if len(data)>=2:
			sp = data[0]
			if sp not in species_kinds:
				species_kinds.append(sp)
				uniq_species.write(line)
		else:
			if sp not in species_kinds:
				species_kinds.append(sp)
				uniq_species.write(line)
print ('ok')
