from pymol import cmd, stored, math
List1=[]

def loadBfacts (mol,startaa=1, visual="Y"):
	obj=cmd.get_object_list(mol)[0]
	cmd.alter(mol,"b=-1.0")
	counter=int(startaa)
	bfacts=[]
	for item in List1:
		bfact=float(item*10)
		bfacts.append(bfact)
		cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
		counter=counter+1
	if visual=="Y":
		cmd.show_as("cartoon",mol)
		# cmd.cartoon("putty", mol)
		cmd.set("cartoon_putty_scale_min", min(bfacts),obj)
		cmd.set("cartoon_putty_scale_max", max(bfacts),obj)
		cmd.set("cartoon_putty_transform", 0,obj)
		cmd.set("cartoon_putty_radius", 0.2,obj)
		cmd.spectrum("b","blue_white_red", "%s and n. CA " %mol)
		cmd.ramp_new("count", obj, [min(bfacts), max(bfacts)],['blue','white','red'])
		cmd.recolor()

cmd.extend("loadBfacts", loadBfacts)
