# This script generates a list of carbohydrate residues based on the ringsize and ring atom names. Then, the coordinates for each ring atom are assigned to variables, and the geometric center is identified for the current residue, and that of the residue connected to C1. Next, the geometric centers of the residues are connected with gray cylinders. Finally, a sphere, cube, cone, or diamond is drawn on top of the current carbohydrate residue. If/then statements are currently used to choose the shape and color for each residue; however, this relies on a list of carbohydrate names which may change over time. Ideally, the sugar identification program developed by Lachele Foley can be used to 1) identify carbohydrate residues, and 2) name those residues solely based on the atomic coordinates of the residue. 

# This script requires the carbohydrate ring to be named as 'C1,C2,C3,C4,C5,O5' - This nomenclature is relatively common, and required by GLYCAM. Note that the ring atoms for sialic acids are C2,C3,C4,C5,C6,O6. Therefore, the current implementation will not work with sialic acid residues.

proc cfg {} {

# This first section can be modified to change the size of the shapes that will be displayed.
	set size 4
	set cube_size [expr $size*0.9]
	set cone_size [expr $size *1.1]
	set star_size $size
	set sphere_size [expr $size/2]
	set diamond_size [expr $size*1.4]

    # Collect a list of residues that contain carbohydrate ring atoms
    set residues [lsort -unique [[atomselect top "ringsize 6 from (hetero and name C1 C2 C3 C4 C5 O5)"] get residue]]
    # For each carbohydrate residue
    foreach residue $residues {

	# IF the residue has atom O1A, it is a sialic acid
	set Sacid [atomselect top "residue $residue and name O1A"]

# For sialic residues, use different atom numbers
        if {[$Sacid get residue] > 0} {
                linked_sialic $residue $diamond_size
        } else {

# For non-sialics, use C1-O5

        # Select the atoms within the carbohydrate ring (note that O5 represents the sixth ring atom) 
        set c1 [atomselect top "residue $residue and name C1"]
        set c2 [atomselect top "residue $residue and name C2"]
        set c3 [atomselect top "residue $residue and name C3"]
        set c4 [atomselect top "residue $residue and name C4"]
        set c5 [atomselect top "residue $residue and name C5"]
        set c6 [atomselect top "residue $residue and name O5"]

        # Get the coordinates for each atom in the 6-membered ring (p stands for point)
        lassign [$c1 get {x y z}] p1
        lassign [$c2 get {x y z}] p2
        lassign [$c3 get {x y z}] p3
        lassign [$c4 get {x y z}] p4
        lassign [$c5 get {x y z}] p5
        lassign [$c6 get {x y z}] p6

	# Calculate the geometric center of the ring
        set ring_atoms [atomselect top "residue $residue and name C1 C2 C3 C4 C5 O5"]
        set geom_center [measure center $ring_atoms]
	# The geometric center of the attached residue is determined later in the script. In some cases, there is no attached residue. The following line assigns the coordinates of the C1 atom to be the 'geometric center of the attached residue' so that the shapes can be aligned properly without crashing when no residue is attached. This occurs for terminal residues not attached to proteins (i.e. OME or OH on reducing terminus).
	set geom_center_att $p1

#############################################################################################################

# Connect the residues with cylinders between geometric centers. This creates a gray cylinder for glycosidic linkage, and a black one to the CA of N-linked or O-linked glycans. This color may change depending on public opinion (option for user?) 

# Identify the geometric center of the ring connected to the C1, or the CA of the amino acid that the glycan is attached to. This geometric center is also used to connect the shapes together via cylinders, and to orient the shape in the correct direction.
                # Variables with _att refer to attributes of the residue that are attached to the C1. Otherwise, the variable refers to the current/targeted residue

        # Identify the atom on the residue connected to C1. O_att and N_att are for attached oxygen and attached nitrogen. N_att is for N-linked glycans
                # O_att stands for attached oxygen. O_att_r is for residue of attached oxygen. O_att_n is for name (i.e. O1).
        set O_att [atomselect top "(oxygen within 2 of residue $residue and name C1) and not name O5"]
                # N_att stands for N-linked 
        set N_att [atomselect top "(nitrogen within 2 of residue $residue and name C1)"]

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
	                        draw color gray
	                        draw cylinder $geom_center $geom_center_att radius 0.5
                        } elseif {[$ring_atoms_att2 get residue] > 0 } {
                                set geom_center_att [measure center $ring_atoms_att2]
                                draw color gray
                                draw cylinder $geom_center $geom_center_att radius 0.5

		# O-linked Glycoproteins: if not a sugar, attached residue is assumed to be a protein. $att_CA stands for Neighboring Protein C-alpha
	                } else {
	                        set att_CA [atomselect top "residue $C_att_r and name CA"]
	                        # If it doesn't have a CA, skip it. This is to avoid crashes on OME residues. CAp stands for CA point
	                        if {[$att_CA get residue] > 0 } {
	                                lassign [$att_CA get {x y z}] att_CAp
	                                draw color black
	                                draw cylinder $geom_center $att_CAp radius 0.5
	                                draw sphere $att_CAp radius 0.5 resolution 36
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
                draw color black
                draw cylinder $geom_center $att_CAp radius 0.5
                draw sphere $att_CAp radius 0.5 resolution 36
			# Note that in this case, $geom_center_att is actually the location of CA. This is only added to align the shapes attached to the protein.
                        lassign [$att_CA get {x y z}] geom_center_att
        }
	
#############################################################################################################

# Depending on the name, change the color/shape (cube, cone, or sphere). This is nowhere near an exhaustive list, and could be edited to follow the logic of the CFG nomenclature (i.e. if the residue has a Nitrogen, draw a cube, etc.)
# It would be nice for the user to have the option to change the size of the shapes. Right now, everything is set to 4 Angstroms, but it's easy to change (check for $size in the other scripts).
set resname [$c1 get resname]

   # The following section augmented by Jodi to contain a more complete list of CFG-supported pyranose sugars (still far from exhaustive)
   # Common || CHARMM36 || GLYCAM06

#source ./ln_oriented_cone.tcl
#source ./ln_oriented_cube.tcl
#source ./ln_oriented_diamond.tcl
#source ./ln_oriented_star.tcl

	# Note:	Currently uses standard VMD colors, such that carbohydrate representations will "match" other default representations.
	# 	Technically, however, CFG nomenclature should use the following colors:
	#	
	#	Yellow (255,255,0)	Currently using "yellow"
	#	Blue (0,0,250)		Currently using "blue"
	#	Green (0,200,50)	Currently using "green"
	#	Red (250,0,0)		Currently using "red"
	#	Orange (250,234,213)	Currently using "yellow"
	#	Purple (200,0,200)
	#	Light blue (233,255,255)
	#	Tan (150,100,50)

	if {$resname == "GAL"  || $resname == "AGAL"  || $resname == "BGAL"  || $resname == "0LA"  || $resname == "0LB"  || $resname == "1LA"  || $resname == "1LB"  || $resname == "2LA"  || $resname == "2LB"  || $resname == "3LA"  || $resname == "3LB"  || $resname == "4LA"  || $resname == "4LB"  || $resname == "6LA"  || $resname == "6LB"  || $resname == "ZLA"  || $resname == "ZLB"  || $resname == "YLA"  || $resname == "YLB"  || $resname == "XLA"  || $resname == "XLB"  || $resname == "WLA"  || $resname == "WLB"  || $resname == "VLA"  || $resname == "VLB"  || $resname == "ULA"  || $resname == "ULB"  || $resname == "TLA"  || $resname == "TLB"  || $resname == "SLA"  || $resname == "SLB"  || $resname == "RLA"  || $resname == "RLB"  || $resname == "QLA"  || $resname == "QLB"  || $resname == "PLA"  || $resname == "PLB" } {
		# Galactose, yellow sphere
		draw color yellow
		draw sphere $geom_center radius $sphere_size resolution 20
	} elseif {$resname == "AGALNA"  || $resname == "BGALNA"  || $resname == "0VA"  || $resname == "0VB"  || $resname == "1VA"  || $resname == "1VB"  || $resname == "2VA"  || $resname == "2VB"  || $resname == "3VA"  || $resname == "3VB"  || $resname == "4VA"  || $resname == "4VB"  || $resname == "6VA"  || $resname == "6VB"  || $resname == "ZVA"  || $resname == "ZVB"  || $resname == "YVA"  || $resname == "YVB"  || $resname == "XVA"  || $resname == "XVB"  || $resname == "WVA"  || $resname == "WVB"  || $resname == "VVA"  || $resname == "VVB"  || $resname == "UVA"  || $resname == "UVB"  || $resname == "TVA"  || $resname == "TVB"  || $resname == "SVA"  || $resname == "SVB"  || $resname == "RVA"  || $resname == "RVB"  || $resname == "QVA"  || $resname == "QVB"  || $resname == "PVA"  || $resname == "PVB" } {
		# N-acetyl-galactosamine, yellow cube
		linked_carb.oriented.cube $p6 $cube_size $geom_center $geom_center_att yellow yellow 
   	} elseif {$resname == "GLC"  || $resname == "AGLC"  || $resname == "BGLC"  || $resname == "0GA"  || $resname == "0GB"  || $resname == "0GA"  || $resname == "0GB"  || $resname == "1GA"  || $resname == "1GB"  || $resname == "2GA"  || $resname == "2GB"  || $resname == "3GA"  || $resname == "3GB"  || $resname == "4GA"  || $resname == "4GB"  || $resname == "6GA"  || $resname == "6GB"  || $resname == "ZGA"  || $resname == "ZGB"  || $resname == "YGA"  || $resname == "YGB"  || $resname == "XGA"  || $resname == "XGB"  || $resname == "WGA"  || $resname == "WGB"  || $resname == "VGA"  || $resname == "VGB"  || $resname == "UGA"  || $resname == "UGB"  || $resname == "TGA"  || $resname == "TGB"  || $resname == "SGA"  || $resname == "SGB"  || $resname == "RGA"  || $resname == "RGB"  || $resname == "QGA"  || $resname == "QGB"  || $resname == "PGA"  || $resname == "PGB" } {
		# Glucose, blue sphere
		draw color blue
		draw sphere $geom_center radius $sphere_size resolution 20
	} elseif {$resname == "NAG"  || $resname == "4YS"  || $resname == "SGN"  || $resname == "AGLCNA"  || $resname == "BGLCNA"  || $resname == "BGLCN0"  || $resname == "0YA"  || $resname == "0YB"  || $resname == "1YA"  || $resname == "1YB"  || $resname == "2YA"  || $resname == "2YB"  || $resname == "3YA"  || $resname == "3YB"  || $resname == "4YA"  || $resname == "4YB"  || $resname == "6YA"  || $resname == "6YB"  || $resname == "ZYA"  || $resname == "ZYB"  || $resname == "YYA"  || $resname == "YYB"  || $resname == "XYA"  || $resname == "XYB"  || $resname == "WYA"  || $resname == "WYB"  || $resname == "VYA"  || $resname == "VYB"  || $resname == "UYA"  || $resname == "UYB"  || $resname == "TYA"  || $resname == "TYB"  || $resname == "SYA"  || $resname == "SYB"  || $resname == "RYA"  || $resname == "RYB"  || $resname == "QYA"  || $resname == "QYB"  || $resname == "PYA"  || $resname == "PYB" } {
		# N-acetyl-glucosamine, blue cube
		linked_carb.oriented.cube $p6 $cube_size $geom_center $geom_center_att blue blue
	} elseif {$resname == "MAN"  || $resname == "BMA"  || $resname == "AMAN"  || $resname == "BMAN"  || $resname == "0MA"  || $resname == "0MB"  || $resname == "1MA"  || $resname == "1MB"  || $resname == "2MA"  || $resname == "2MB"  || $resname == "3MA"  || $resname == "3MB"  || $resname == "4MA"  || $resname == "4MB"  || $resname == "6MA"  || $resname == "6MB"  || $resname == "ZMA"  || $resname == "ZMB"  || $resname == "YMA"  || $resname == "YMB"  || $resname == "XMA"  || $resname == "XMB"  || $resname == "WMA"  || $resname == "WMB"  || $resname == "VMA"  || $resname == "VMB"  || $resname == "UMA"  || $resname == "UMB"  || $resname == "TMA"  || $resname == "TMB"  || $resname == "SMA"  || $resname == "SMB"  || $resname == "RMA"  || $resname == "RMB"  || $resname == "QMA"  || $resname == "QMB"  || $resname == "PMA"  || $resname == "PMB" } {
		# Mannose, green sphere 
		draw color green
		draw sphere $geom_center radius $sphere_size resolution 20
	} elseif {$resname == "0WA"  || $resname == "0WB"  || $resname == "1WA"  || $resname == "1WB"  || $resname == "2WA"  || $resname == "2WB"  || $resname == "3WA"  || $resname == "3WB"  || $resname == "4WA"  || $resname == "4WB"  || $resname == "6WA"  || $resname == "6WB"  || $resname == "ZWA"  || $resname == "ZWB"  || $resname == "YWA"  || $resname == "YWB"  || $resname == "XWA"  || $resname == "XWB"  || $resname == "WWA"  || $resname == "WWB"  || $resname == "VWA"  || $resname == "VWB"  || $resname == "UWA"  || $resname == "UWB"  || $resname == "TWA"  || $resname == "TWB"  || $resname == "SWA"  || $resname == "SWB"  || $resname == "RWA"  || $resname == "RWB"  || $resname == "QWA"  || $resname == "QWB"  || $resname == "PWA"  || $resname == "PWB" } {
		# N-acetyl-mannosamine, green cube 
		linked_carb.oriented.cube $p6 $cube_size $geom_center $geom_center_att green green
	} elseif {$resname == "FUC"  || $resname == "AFUC"  || $resname == "BFUC"  || $resname == "0fA"  || $resname == "0fB"  || $resname == "1fA"  || $resname == "1fB"  || $resname == "2fA"  || $resname == "2fB"  || $resname == "3fA"  || $resname == "3fB"  || $resname == "4fA"  || $resname == "4fB"  || $resname == "6fA"  || $resname == "6fB"  || $resname == "ZfA"  || $resname == "ZfB"  || $resname == "YfA"  || $resname == "YfB"  || $resname == "XfA"  || $resname == "XfB"  || $resname == "WfA"  || $resname == "WfB"  || $resname == "VfA"  || $resname == "VfB"  || $resname == "UfA"  || $resname == "UfB"  || $resname == "TfA"  || $resname == "TfB"  || $resname == "SfA"  || $resname == "SfB"  || $resname == "RfA"  || $resname == "RfB"  || $resname == "QfA"  || $resname == "QfB"  || $resname == "PfA"  || $resname == "PfB" } {
		# Fucose, red cone
		draw color red
		linked_carb.oriented.cone $p6 $cone_size $geom_center $geom_center_att red red
	} elseif {$resname == "XYS"  || $resname == "LXC"  || $resname == "AXYL"  || $resname == "BXYL"  || $resname == "0XA"  || $resname == "0XB"  || $resname == "1XA"  || $resname == "1XB"  || $resname == "2XA"  || $resname == "2XB"  || $resname == "3XA"  || $resname == "3XB"  || $resname == "4XA"  || $resname == "4XB"  || $resname == "6XA"  || $resname == "6XB"  || $resname == "ZXA"  || $resname == "ZXB"  || $resname == "YXA"  || $resname == "YXB"  || $resname == "XXA"  || $resname == "XXB"  || $resname == "WXA"  || $resname == "WXB"  || $resname == "VXA"  || $resname == "VXB"  || $resname == "UXA"  || $resname == "UXB"  || $resname == "TXA"  || $resname == "TXB"  || $resname == "SXA"  || $resname == "SXB"  || $resname == "RXA"  || $resname == "RXB"  || $resname == "QXA"  || $resname == "QXB"  || $resname == "PXA"  || $resname == "PXB" } {
		# Xylose, orange star
		draw color orange2
                linked_carb.oriented.star $p6 $star_size $geom_center $geom_center_att
	} elseif {$resname == "AGLCA"  || $resname == "BGLCA"  || $resname == "BGLCA0"  || $resname == "0ZA"  || $resname == "0ZB"  || $resname == "1ZA"  || $resname == "1ZB"  || $resname == "2ZA"  || $resname == "2ZB"  || $resname == "3ZA"  || $resname == "3ZB"  || $resname == "4ZA"  || $resname == "4ZB"  || $resname == "6ZA"  || $resname == "6ZB"  || $resname == "ZZA"  || $resname == "ZZB"  || $resname == "YZA"  || $resname == "YZB"  || $resname == "XZA"  || $resname == "XZB"  || $resname == "WZA"  || $resname == "WZB"  || $resname == "VZA"  || $resname == "VZB"  || $resname == "UZA"  || $resname == "UZB"  || $resname == "TZA"  || $resname == "TZB"  || $resname == "SZA"  || $resname == "SZB"  || $resname == "RZA"  || $resname == "RZB"  || $resname == "QZA"  || $resname == "QZB"  || $resname == "PZA"  || $resname == "PZB" } {
		# Glucuronic acid, blue/white diamond
		linked_carb.oriented.diamond $p6 $diamond_size $geom_center $geom_center_att blue white
	} elseif {$resname == "IDS"  || $resname == "AIDOA"  || $resname == "BIDOA"  || $resname == "0uA"  || $resname == "0uB"  || $resname == "1uA"  || $resname == "1uB"  || $resname == "2uA"  || $resname == "2uB"  || $resname == "3uA"  || $resname == "3uB"  || $resname == "4uA"  || $resname == "4uB"  || $resname == "6uA"  || $resname == "6uB"  || $resname == "ZuA"  || $resname == "ZuB"  || $resname == "YuA"  || $resname == "YuB"  || $resname == "XuA"  || $resname == "XuB"  || $resname == "WuA"  || $resname == "WuB"  || $resname == "VuA"  || $resname == "VuB"  || $resname == "UuA"  || $resname == "UuB"  || $resname == "TuA"  || $resname == "TuB"  || $resname == "SuA"  || $resname == "SuB"  || $resname == "RuA"  || $resname == "RuB"  || $resname == "QuA"  || $resname == "QuB"  || $resname == "PuA"  || $resname == "PuB" } {
		# L-Iduronic acid, white/tan diamond
		linked_carb.oriented.diamond $p6 $diamond_size $geom_center $geom_center_att white ochre
	} elseif {$resname == "0OA"  || $resname == "0OB"  || $resname == "1OA"  || $resname == "1OB"  || $resname == "2OA"  || $resname == "2OB"  || $resname == "3OA"  || $resname == "3OB"  || $resname == "4OA"  || $resname == "4OB"  || $resname == "6OA"  || $resname == "6OB"  || $resname == "ZOA"  || $resname == "ZOB"  || $resname == "YOA"  || $resname == "YOB"  || $resname == "XOA"  || $resname == "XOB"  || $resname == "WOA"  || $resname == "WOB"  || $resname == "VOA"  || $resname == "VOB"  || $resname == "UOA"  || $resname == "UOB"  || $resname == "TOA"  || $resname == "TOB"  || $resname == "SOA"  || $resname == "SOB"  || $resname == "ROA"  || $resname == "ROB"  || $resname == "QOA"  || $resname == "QOB"  || $resname == "POA"  || $resname == "POB" } {
		# GalA, yellow|white diamond
		linked_carb.oriented.diamond $p6 $diamond_size $geom_center $geom_center_att yellow white
        } elseif {$resname == "RAM"  || $resname == "ARHM"  || $resname == "BRHM"  || $resname == "0HA"  || $resname == "0HB"  || $resname == "1HA"  || $resname == "1HB"  || $resname == "2HA"  || $resname == "2HB"  || $resname == "3HA"  || $resname == "3HB"  || $resname == "4HA"  || $resname == "4HB"  || $resname == "6HA"  || $resname == "6HB"  || $resname == "ZHA"  || $resname == "ZHB"  || $resname == "YHA"  || $resname == "YHB"  || $resname == "XHA"  || $resname == "XHB"  || $resname == "WHA"  || $resname == "WHB"  || $resname == "VHA"  || $resname == "VHB"  || $resname == "UHA"  || $resname == "UHB"  || $resname == "THA"  || $resname == "THB"  || $resname == "SHA"  || $resname == "SHB"  || $resname == "RHA"  || $resname == "RHB"  || $resname == "QHA"  || $resname == "QHB"  || $resname == "PHA"  || $resname == "PHB" } {
		# Rhamnose, gray cone
		draw color gray
		linked_carb.oriented.cone $p6 $cone_size $geom_center $geom_center_att gray gray
	} elseif {$resname == "UYS" || $resname == "WYS"} {
		# Some GAG support: Sulfated glucosamine, white\blue cube
		linked_carb.oriented.cube $p6 $cube_size $geom_center $geom_center_att white blue
	} else {
		# Color by Resname if not listed above, allowing the user the option to change the color in the GUI
	        draw color [colorinfo category Resname $resname]
                draw sphere $geom_center radius $sphere_size resolution 20
	}
	}
    }
}
