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

proc linked_carb.oriented.cube {p6 size geom_center geom_center_att color1 color2} {

        # The objective is to orient the shapes so that they face the neighboring residue, connected at C1. This is accomplished by calculating the equation for the line that connects the geometric centers of the target residue (where the shape will be placed) and the neighboring residue (which is attached at the C1 atom). If we consider the geometric center of the target residue to be point A, the center of the neighboring residue to be point B, and the origin to be point O, then vector_AB=vector_OB - vector_OA. 
        # Knowing vector AB allows us to put new points along the line; however, we want the shapes to be of a defined size. In order to adjust the size properly, the distance between geometric centers is determined, scaled to match the desired size, and then added/subtracted from the geometric center of the target sugar.

        set vec_AB [vecsub $geom_center_att $geom_center]
                # Resize the shape. $size refers to the total size, $half_length refers to the distance required create points on either side of the geometric center.
        set half_length [expr $size/2]
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

####################### This part is different from drawing a cone #######################

		# perp_1 represents a point perpendicular to the previously created two
	set perp_1 [vecnorm [veccross [vecsub $x1 $x2] [vecsub $x2 $p6]]]
        set perp1 [vecscale $perp_1 $half_length]
        set perp2 [vecscale $perp_1 -$half_length]
	        # Each 'o' represents a point on the box (o stands for original points, which are based on the two that lie along the line connecting the two residues)
        set o1 [vecadd $perp1 $x1]
        set o2 [vecadd $perp1 $x2]
        set o3 [vecadd $perp2 $x1]
        set o4 [vecadd $perp2 $x2]

# This is only used for testing, to see where the points are
#	draw color black
#	draw sphere $o1 radius 0.2
#	draw sphere $o2 radius 0.2
#	draw sphere $o3 radius 0.2
#	draw sphere $o4 radius 0.2

                # This creates 8 points (the corners of the box) based upon the coordinates of o1-o4
        set perp_for [vecscale [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]] $half_length]
        set perp_back [vecscale [vecnorm [veccross [vecsub $o1 $o2] [vecsub $o3 $o1]]] -$half_length]
        set s1 [vecadd $perp_for $o1]
        set s2 [vecadd $perp_for $o2]
        set s3 [vecadd $perp_for $o3]
        set s4 [vecadd $perp_for $o4]
        set s5 [vecadd $perp_back $o1]
        set s6 [vecadd $perp_back $o2]
        set s7 [vecadd $perp_back $o3]
        set s8 [vecadd $perp_back $o4]

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

        # Draw the Cube - some cubes require two colors, so the triangles are intentionally divided so that one side shows both colors
        # Color 1 (blue, yellow, or green)
        draw color $color1
        draw triangle $s2 $s3 $s4
        draw triangle $s1 $s2 $s6
        draw triangle $s4 $s8 $s7
        draw triangle $s5 $s6 $s8
        draw triangle $s2 $s4 $s8
        draw triangle $s1 $s5 $s7
        # Color 2 (white)
        draw color $color2
        draw triangle $s1 $s2 $s3
        draw triangle $s3 $s4 $s7
        draw triangle $s1 $s5 $s6
        draw triangle $s5 $s7 $s8
        draw triangle $s2 $s6 $s8
        draw triangle $s1 $s3 $s7

        # Draw a border around the edges of the cube - this is important for white shapes that blend into white backgrounds, and should be an option for the user
#       draw color gray
#       draw cylinder $s2 $s4 radius 0.01
#       draw cylinder $s4 $s8 radius 0.01
#       draw cylinder $s8 $s6 radius 0.01
#       draw cylinder $s6 $s2 radius 0.01
#       draw cylinder $s1 $s3 radius 0.01
#       draw cylinder $s3 $s7 radius 0.01
#       draw cylinder $s7 $s5 radius 0.01
#       draw cylinder $s5 $s1 radius 0.01
#       draw cylinder $s3 $s4 radius 0.01
#       draw cylinder $s8 $s7 radius 0.01
#       draw cylinder $s6 $s5 radius 0.01
#       draw cylinder $s1 $s2 radius 0.01

}
