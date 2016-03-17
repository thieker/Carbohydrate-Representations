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

proc linked_carb.oriented.diamond {p6 size geom_center geom_center_att color1 color2} {

	# The objective is to orient the shapes so that they face the neighboring residue, connected at C1. This is accomplished by calculating the equation for the line that connects the geometric centers of the target residue (where the shape will be placed) and the neighboring residue (which is attached at the C1 atom). If we consider the geometric center of the target residue to be point A, the center of the neighboring residue to be point B, and the origin to be point O, then vector_AB=vector_OB - vector_OA. 
	# Knowing vector AB allows us to put new points along the line; however, we want the shapes to be of a defined size. In order to adjust the size properly, the distance between geometric centers is determined, scaled to match the desired size, and then added/subtracted from the geometric center of the target sugar.

        set vec_AB [vecsub $geom_center_att $geom_center]
                # Resize the shape. $size refers to the total size, $diamond_size refers to the distance required create points on either side of the geometric center.
        set diamond_size [expr $size*0.5]
		# Vec_AB is currently too large, we want it to be adjusted so that the distance is equal to the offset that will be used to place the ponits around the geometric center. The equation asks the question, what factor should be multiplied times the distance in order to prodce $half-length?
	set adjustment [expr $diamond_size/[vecdist $geom_center $geom_center_att]]
		# Adjust vector_AB by the amount determined in both the forward and reverse directions
	set adj_vec_AB_1 [vecscale $adjustment $vec_AB]
	set adj_vec_AB_2 [vecscale -$adjustment $vec_AB]
		# Add two points along the line connecting the residues in the forward and reverse direction
	set x1 [vec + $adj_vec_AB_1 $geom_center]
	set x2 [vec + $adj_vec_AB_2 $geom_center]

                # perp_1 represents a point perpendicular to the previously created two
        set perp_1 [vecnorm [veccross [vecsub $x1 $x2] [vecsub $x2 $p6]]]
        set perp1 [vecscale $perp_1 $diamond_size]
        set perp2 [vecscale $perp_1 -$diamond_size]

		# The number of sides to the diamond is up for debate. It was suggested to use six sides since the diamond could look like a cube that has been rotated; however, four seems to work since the cylinder goes directly through one corner. Maybe it could be an option for the user?
                # Each 'o' represents a corner of a square that is centered on $geom_center
        set o1 [vecadd $perp1 $x1]
        set o2 [vecadd $perp1 $x2]
        set o3 [vecadd $perp2 $x1]
        set o4 [vecadd $perp2 $x2]

        # Determine coordinates for outer points of the diamond
                # Top point of star is $diamond_size above geometric center
                # Vector distance between center and top point of star
        set d1 [vecsub $x2 $geom_center]
                # Rotate by 72 degrees to identify other points of the star - note that t5 is the same as n1
        set t1 [vectrans [trans angle $o1 $geom_center $o2 90] $d1]
        set t2 [vectrans [trans angle $o1 $geom_center $o2 180] $d1]
        set t3 [vectrans [trans angle $o1 $geom_center $o2 270] $d1]
        set t4 [vectrans [trans angle $o1 $geom_center $o2 360] $d1]
        set outer_1 [vecadd $t1 $geom_center]
        set outer_2 [vecadd $t2 $geom_center]
        set outer_3 [vecadd $t3 $geom_center]
        set outer_4 [vecadd $t4 $geom_center]

               # The following function creates the top and bottom of the diamond by creating two points at the geometric center that are perpendicular to the plane of the square (and parallel with the plane of the ring).
        set perp_for [vecscale $diamond_size [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]]]
        set perp_back [vecscale -$diamond_size [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]]]
        set top [vecadd $perp_for $geom_center]
        set bottom [vecadd $perp_back $geom_center]

# Draw the diamond
        # Front
	draw color $color1
        draw triangle $outer_1 $outer_2 $bottom
        draw triangle $outer_1 $outer_4 $bottom
        draw triangle $outer_3 $outer_2 $top
        draw triangle $outer_3 $outer_4 $top
        # Back
	draw color $color2
        draw triangle $outer_1 $outer_2 $top
        draw triangle $outer_1 $outer_4 $top
        draw triangle $outer_3 $outer_2 $bottom
        draw triangle $outer_3 $outer_4 $bottom
        # Connect
#        draw triangle $outer_1 $outer_4 $outer_1


#	draw color $color1
#        draw cone $geom_center $x1 radius $diamond_size resolution 4
#	draw color $color2
#        draw cone $geom_center $x2 radius $diamond_size resolution 4
		# Only used to ensure that the shape is centered properly, and of the appropriate size
#	draw sphere $geom_center radius 1.75 resolution 40
}
