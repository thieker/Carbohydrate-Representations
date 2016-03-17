########################################################################################################################################################

# This section allows for vector math functions, and was obtained from http://wiki.tcl.tk/14022
 proc lmap {_var list body} {
     upvar 1 $_var var
     set res {}
     foreach var $list {lappend res [uplevel 1 $body]}
     set res
 }
#-- We need basic scalar operators from [expr] factored out:
 foreach op {+ - * / % ==} {proc $op {a b} "expr {\$a $op \$b}"}

 proc vec {op a b} {
     if {[llength $a] == 1 && [llength $b] == 1} {
         $op $a $b
     } elseif {[llength $a]==1} {
         lmap i $b {vec $op $a $i}
     } elseif {[llength $b]==1} {
         lmap i $a {vec $op $i $b}
     } elseif {[llength $a] == [llength $b]} {
         set res {}
         foreach i $a j $b {lappend res [vec $op $i $j]}
         set res
     } else {error "length mismatch [llength $a] != [llength $b]"}
 }

########################################################################################################################################################

proc linked_carb.oriented.star {p6 size geom_center geom_center_att} {

        # The objective is to orient the shapes so that they face the neighboring residue, connected at C1. This is accomplished by calculating the equation for the line that connects the geometric centers of the target residue (where the shape will be placed) and the neighboring residue (which is attached at the C1 atom). If we consider the geometric center of the target residue to be point A, the center of the neighboring residue to be point B, and the origin to be point O, then vector_AB=vector_OB - vector_OA. 
        # Knowing vector AB allows us to put new points along the line; however, we want the shapes to be of a defined size. In order to adjust the size properly, the distance between geometric centers is determined, scaled to match the desired size, and then added/subtracted from the geometric center of the target sugar.

        set vec_AB [vecsub $geom_center_att $geom_center]
                # Resize the shape. $size refers to the total size, $half_length refers to the distance required create points on either side of the geometric center.
	set star_size [expr $size*1.5]
        set half_length [expr $star_size/2]
        set thickness [expr $star_size/4]
                # Vec_AB is currently too large, we want it to be adjusted so that the distance is equal to the offset that will be used to place the ponits around the geometric center. The equation asks the question, what factor should be multiplied times the distance in order to prodce $half-length?
        set adjustment [expr $half_length/[vecdist $geom_center $geom_center_att]]
                # Adjust vector_AB by the amount determined in both the forward and reverse directions
        set adj_vec_AB_1 [vecscale $adjustment $vec_AB]
        set adj_vec_AB_2 [vecscale -$adjustment $vec_AB]
                # Add two points along the line connecting the residues in the forward and reverse direction
        set x1 [vec + $adj_vec_AB_1 $geom_center]
        set x2 [vec + $adj_vec_AB_2 $geom_center]

# This is only used for testing, to see where the points are
#        draw color red
#	draw sphere $x1 radius 0.2
#	draw sphere $x2 radius 0.2
		# perp_1 represents a point perpendicular to the previously created two
	set perp_1 [vecnorm [veccross [vecsub $x1 $x2] [vecsub $x2 $p6]]]
        set perp1 [vecscale $perp_1 $half_length]
        set perp2 [vecscale $perp_1 -$half_length]
	        # Each 'o' represents a corner of a square that is centered on $geom_center
        set o1 [vecadd $perp1 $x1]
        set o2 [vecadd $perp1 $x2]
        set o3 [vecadd $perp2 $x1]
        set o4 [vecadd $perp2 $x2]

	# Determine coordinates for outer points of the star
		# Top point of star is $half_length above geometric center
		# Vector distance between center and top point of star
	set d1 [vecsub $x2 $geom_center]
		# Rotate by 72 degrees to identify other points of the star 
	set t1 [vectrans [trans angle $o1 $geom_center $o2 72] $d1]
	set t2 [vectrans [trans angle $o1 $geom_center $o2 144] $d1]
	set t3 [vectrans [trans angle $o1 $geom_center $o2 216] $d1]
	set t4 [vectrans [trans angle $o1 $geom_center $o2 288] $d1]
	set t5 [vectrans [trans angle $o1 $geom_center $o2 360] $d1]
	set outer_1 [vecadd $t1 $geom_center]
	set outer_2 [vecadd $t2 $geom_center]
	set outer_3 [vecadd $t3 $geom_center]
	set outer_4 [vecadd $t4 $geom_center]
	set outer_5 [vecadd $t5 $geom_center]

	# Determine coordinates for inner points of the star
		# Bottom point of star is $half_length below geometric center
		# Vector distance between center and top point of star
	set d2 [vecscale [vecsub $x1 $geom_center] 0.5]
		# Rotate by 72 degrees to identify other points of the star - note that b5 is the same as n2
	set b1 [vectrans [trans angle $o1 $geom_center $o2 72] $d2]
	set b2 [vectrans [trans angle $o1 $geom_center $o2 144] $d2]
	set b3 [vectrans [trans angle $o1 $geom_center $o2 216] $d2]
	set b4 [vectrans [trans angle $o1 $geom_center $o2 288] $d2]
	set b5 [vectrans [trans angle $o1 $geom_center $o2 360] $d2]
	set inner_1 [vecadd $b1 $geom_center]
	set inner_2 [vecadd $b2 $geom_center]
	set inner_3 [vecadd $b3 $geom_center]
	set inner_4 [vecadd $b4 $geom_center]
	set inner_5 [vecadd $b5 $geom_center]
	
                # The following function makes the star 3D by creating two points at the geometric center that are perpendicular to the plane of the star
        set perp_for [vecscale $thickness [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]]]
        set perp_back [vecscale -$thickness [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]]]
        set center_1 [vecadd $perp_for $geom_center]
        set center_2 [vecadd $perp_back $geom_center]

# Draw the star
	draw triangle $outer_1 $inner_3 $center_1
	draw triangle $outer_1 $inner_3 $center_2
	draw triangle $outer_1 $inner_4 $center_1
	draw triangle $outer_1 $inner_4 $center_2
	draw triangle $outer_2 $inner_4 $center_1
	draw triangle $outer_2 $inner_4 $center_2
	draw triangle $outer_2 $inner_5 $center_1
	draw triangle $outer_2 $inner_5 $center_2
	draw triangle $outer_3 $inner_5 $center_1
	draw triangle $outer_3 $inner_5 $center_2
	draw triangle $outer_3 $inner_1 $center_1
	draw triangle $outer_3 $inner_1 $center_2
	draw triangle $outer_4 $inner_1 $center_1
	draw triangle $outer_4 $inner_1 $center_2
	draw triangle $outer_4 $inner_2 $center_1
	draw triangle $outer_4 $inner_2 $center_2
	draw triangle $outer_5 $inner_2 $center_1
	draw triangle $outer_5 $inner_2 $center_2
	draw triangle $outer_5 $inner_3 $center_1
	draw triangle $outer_5 $inner_3 $center_2

######################################################################################################3

# This is only used for testing, to see where the points are
#	draw color black
#	draw sphere $o1 radius 0.2
#	draw color yellow
#	draw sphere $o2 radius 0.2
#	draw color green
#	draw sphere $o3 radius 0.2
#	draw color orange
#	draw sphere $o4 radius 0.2
#	draw color purple
#	draw sphere $n1 radius 0.2
#	draw color green
#	draw sphere $outer_1 radius 0.2
#	draw sphere $corner_1 radius 0.2
#	draw color blue
#	draw sphere $outer_2 radius 0.2
#	draw sphere $corner_2 radius 0.2
#	draw color yellow
#	draw sphere $outer_3 radius 0.2
#	draw color ochre
#	draw sphere $outer_4 radius 0.2
#	draw color orange
#	draw sphere $outer_5 radius 0.2
#	draw color green
#	draw sphere $inner_1 radius 0.2
#	draw color blue
#	draw sphere $inner_2 radius 0.2
#	draw color yellow
#	draw sphere $inner_3 radius 0.2
#	draw color ochre
#	draw sphere $inner_4 radius 0.2
#	draw color orange
#	draw sphere $inner_5 radius 0.2
	
# This is only used for testing, to see where the points are
#	draw color yellow
#	draw sphere $s1 radius 0.2
#	draw sphere $s2 radius 0.2
#	draw sphere $s3 radius 0.2
#	draw sphere $s4 radius 0.2
#	draw sphere $s5 radius 0.2
#	draw sphere $s6 radius 0.2
#	draw sphere $s7 radius 0.2
#	draw sphere $s8 radius 0.2
#	draw color blue

}
