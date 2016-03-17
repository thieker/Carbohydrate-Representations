# This script generates a list of carbohydrate residues based on the ringsize and ring atom names. Then, the coordinates for each ring atom are assigned to variables, and the geometric center is identified for the current residue, and that of the residue connected to C1. Next, the geometric centers of the residues are connected with gray cylinders. Finally, a sphere, cube, cone, or diamond is drawn on top of the current carbohydrate residue. If/then statements are currently used to choose the shape and color for each residue; however, this relies on a list of carbohydrate names which may change over time. Ideally, the sugar identification program developed by Lachele Foley can be used to 1) identify carbohydrate residues, and 2) name those residues solely based on the atomic coordinates of the residue. 

# This script requires the carbohydrate ring to be named as 'C1,C2,C3,C4,C5,O5' - This nomenclature is relatively common, and required by GLYCAM. Note that the ring atoms for sialic acids are C2,C3,C4,C5,C6,O6. Therefore, the current implementation will not work with sialic acid residues.

proc linked_sialic {residue size} {

        # Select the atoms within the carbohydrate ring (note c6 represents the O5 atom)
        set c1 [atomselect top "residue $residue and name C2"]
        set c2 [atomselect top "residue $residue and name C3"]
        set c3 [atomselect top "residue $residue and name C4"]
        set c4 [atomselect top "residue $residue and name C5"]
        set c5 [atomselect top "residue $residue and name C6"]
        set c6 [atomselect top "residue $residue and name O6"]

        # Get the coordinates for each atom in the 6-membered ring (p stands for point)
        lassign [$c1 get {x y z}] p1
        lassign [$c2 get {x y z}] p2
        lassign [$c3 get {x y z}] p3
        lassign [$c4 get {x y z}] p4
        lassign [$c5 get {x y z}] p5
        lassign [$c6 get {x y z}] p6

	# Calculate the geometric center of the ring
        set ring_atoms [atomselect top "residue $residue and name C2 C3 C4 C5 C6 O6"]
        set geom_center [measure center $ring_atoms]
	# The geometric center of the attached residue is determined later in the script. In some cases, there is no attached residue. The following line assigns the coordinates of the C1 atom to be the 'geometric center of the attached residue' so that the shapes can be aligned properly without crashing when no residue is attached. This occurs for terminal residues not attached to proteins (i.e. OME or OH on reducing terminus).
	set geom_center_att $p1
#############################################################################################################

# Connect the residues with cylinders between geometric centers. This creates a gray cylinder for glycosidic linkage, and a black one to the CA of N-linked or O-linked glycans. This color may change depending on public opinion (option for user?) 

# Identify the geometric center of the ring connected to the C1, or the CA of the amino acid that the glycan is attached to. This geometric center is also used to connect the shapes together via cylinders, and to orient the shape in the correct direction.
                # Variables with _att refer to attributes of the residue that are attached to the C1. Otherwise, the variable refers to the current/targeted residue

        # Identify the atom on the residue connected to C1. O_att and N_att are for attached oxygen and attached nitrogen. N_att is for N-linked glycans
                # O_att stands for attached oxygen. O_att_r is for residue of attached oxygen. O_att_n is for name (i.e. O1). Note that the atom selection is different than non-sialic acids because O6 is a possible hydroxyl group for attached monosaccharides.
        set O_att [atomselect top "(oxygen within 2 of residue $residue and name C2) and not residue $residue"]
                # N_att stands for N-linked 
        set N_att [atomselect top "(nitrogen within 2 of residue $residue and name C2)"]

        # If the residue has an attached oxygen, draw a cylinder. 
        if {[$O_att get residue] > 0 } {
                        #O_att_r is for residue of attached oxygen. O_att_n is for atom name (i.e. O1).
                lassign [$O_att get residue] O_att_r
                lassign [$O_att get name] O_att_n
                        # C_att stands for carbon of attached residue
                set C_att [atomselect top "(carbon within 2 of residue $O_att_r and name $O_att_n) and not residue $residue"]
		if {[$C_att get residue] > 0 } {
	                lassign [$C_att get residue] C_att_r

		# Glycosidic Linkage (attached residue is a sugar)
	                        # ring_atoms_att refers to the ring atoms of the attached carbohydrate residue, ring_atoms_att2 is for sialic acids
	                set ring_atoms_att [atomselect top "residue $C_att_r and name C1 C2 C3 C4 C5 O5"]
	                set ring_atoms_att2 [atomselect top "residue $C_att_r and name C2 C3 C4 C5 C6 O6"]
	                if {[$ring_atoms_att get residue] > 0 } {
	                        set geom_center_att [measure center $ring_atoms_att]
				# If the user desires a small shape, the connecting cylinders will not be displayed
				if {$size > 3 } {
		                        draw color gray
		                        draw cylinder $geom_center $geom_center_att radius 0.5
				} else {}
                        } elseif {[$ring_atoms_att2 get residue] > 0 } {
	                        set geom_center_att [measure center $ring_atoms_att2]
                                if {$size > 3 } {
		                        draw color gray
		                        draw cylinder $geom_center $geom_center_att radius 0.5
				} else {}
		# O-linked Glycoproteins: if not a sugar, attached residue is assumed to be a protein. $att_CA stands for Neighboring Protein C-alpha
	                } else {
	                        set att_CA [atomselect top "residue $C_att_r and name CA"]
	                        # If it doesn't have a CA, skip it. This is to avoid crashes on OME residues. CAp stands for CA point
	                        if {[$att_CA get residue] > 0 } {
	                                lassign [$att_CA get {x y z}] att_CAp
	                                if {$size > 4 } {
		                                draw color black
		                                draw cylinder $geom_center $att_CAp radius 0.5
		                                draw sphere $att_CAp radius 0.5 resolution 36
					} else {}
					# Note that in this case, $geom_center_att is actually the location of CA. This is necessary to align the shapes later
	                                lassign [$att_CA get {x y z}] geom_center_att
				}
	                }
		} 
	}
	
	# N-linked Glycoproteins: if a nitrogen is within 2 of C1, the attached residue is assumed to be a protein
        if {[$N_att get residue] > 0 } {
                lassign [$N_att get residue] N_att_r
                set att_CA [atomselect top "residue $N_att_r and name CA"]
                lassign [$att_CA get {x y z}] att_CAp
                if {$size > 4 } {
	                draw color black
	                draw cylinder $geom_center $att_CAp radius 0.5
	                draw sphere $att_CAp radius 0.5 resolution 36
		} else {}
			# Note that in this case, $geom_center_att is actually the location of CA. This is only added to align the shapes attached to the protein.
                        lassign [$att_CA get {x y z}] geom_center_att
        }
	
#############################################################################################################

# Depending on the name, change the color/shape (cube, cone, or sphere). This is nowhere near an exhaustive list, and could be edited to follow the logic of the CFG nomenclature (i.e. if the residue has a Nitrogen, draw a cube, etc.)
# It would be nice for the user to have the option to change the size of the shapes. Right now, everything is set to 4 Angstroms, but it's easy to change (check for $size in the other scripts).

		linked_carb.oriented.diamond $p6 $size $geom_center $geom_center_att violet2 violet2
}
