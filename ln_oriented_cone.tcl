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

proc linked_carb.oriented.cone {p6 size geom_center geom_center_att color1 color2} {

        # The objective is to orient the shapes so that they face the neighboring residue, connected at C1. This is accomplished by calculating the equation for the line that connects the geometric centers of the target residue (where the shape will be placed) and the neighboring residue (which is attached at the C1 atom). If we consider the geometric center of the target residue to be point A, the center of the neighboring residue to be point B, and the origin to be point O, then vector_AB=vector_OB - vector_OA. 
        # Knowing vector AB allows us to put new points along the line; however, we want the shapes to be of a defined size. In order to adjust the size properly, the distance between geometric centers is determined, scaled to match the desired size, and then added/subtracted from the geometric center of the target sugar.

        set vec_AB [vecsub $geom_center_att $geom_center]
                # Resize the shape. $size refers to the total size, $half_length refers to the distance required create points on either side of the geometric center.
        set half_length [expr $size/2]
                # Vec_AB is currently too large, we want it to be adjusted so that the distance is equal to the offset that will be used to place the points around the geometric center. The equation asks the question, what factor should be multiplied times the distance in order to prodce $half-length?
                # Two adjustments are added to shift the geom_center of the shape.
        set adjustment1 [expr ($half_length*0.66)/[vecdist $geom_center $geom_center_att]]
        set adjustment2 [expr ($half_length*1.33)/[vecdist $geom_center $geom_center_att]]
                # Adjust vector_AB by the amount determined in both the forward and reverse directions
#       set adj_vec_AB_1 [vecscale $adjustment $vec_AB]
        set adj_vec_AB_1 [vecscale $adjustment1 $vec_AB]
        set adj_vec_AB_2 [vecscale -$adjustment2 $vec_AB]
                # Add two points along the line connecting the residues in the forward and reverse direction
        set x1 [vecadd $adj_vec_AB_1 $geom_center]
        set x2 [vecadd $adj_vec_AB_2 $geom_center]

# This is only used for testing, to see where the points are
#        draw color red
#       draw sphere $x1 radius 0.2
#       draw sphere $x2 radius 0.2
                # perp_1 represents a point perpendicular to the previously created two
        set perp_1 [vecnorm [veccross [vecsub $x1 $x2] [vecsub $x2 $geom_center_att]]]
        set perp1 [vecscale $perp_1 1]
        set perp2 [vecscale $perp_1 -0.5]
        set o1 [vecadd $perp1 $x1]
        set o2 [vecadd $perp2 $x1]

        set perp_3 [vecnorm [veccross [vecsub $o1 $x2] [vecsub $x2 $geom_center_att]]]
        set perp3 [vecscale $perp_3 $half_length]
        set o3 [vecadd $perp3 $x1]

        # Determine coordinates for outer points of the rectangle
                # Top point of rectangle is $half_length above geometric center
                # Vector distance between center and top point of rectangle
        set d1 [vecsub $o3 $x1]
                # Rotate by 90 degrees to identify other points of the rectangle - note that t5 is the same as x1
        set t1 [vectrans [trans angle $o1 $o3 $x1 45] $d1]
        set t2 [vectrans [trans angle $o1 $o3 $x1 90] $d1]
        set t3 [vectrans [trans angle $o1 $o3 $x1 135] $d1]
        set t4 [vectrans [trans angle $o1 $o3 $x1 180] $d1]
        set t5 [vectrans [trans angle $o1 $o3 $x1 225] $d1]
        set t6 [vectrans [trans angle $o1 $o3 $x1 270] $d1]
        set t7 [vectrans [trans angle $o1 $o3 $x1 315] $d1]
        set t8 [vectrans [trans angle $o1 $o3 $x1 360] $d1]
        set outer_1 [vecadd $t1 $x1]
        set outer_2 [vecadd $t2 $x1]
        set outer_3 [vecadd $t3 $x1]
        set outer_4 [vecadd $t4 $x1]
        set outer_5 [vecadd $t5 $x1]
        set outer_6 [vecadd $t6 $x1]
        set outer_7 [vecadd $t7 $x1]
        set outer_8 [vecadd $t8 $x1]

# Draw the rectangle
        # Base
	draw color $color1
        draw triangle $outer_1 $outer_2 $x1
        draw triangle $outer_1 $outer_8 $x1
        draw triangle $outer_3 $outer_2 $x1
        draw triangle $outer_3 $outer_4 $x1
	draw color $color2
        draw triangle $outer_5 $outer_4 $x1
        draw triangle $outer_5 $outer_6 $x1
        draw triangle $outer_7 $outer_6 $x1
        draw triangle $outer_7 $outer_8 $x1

	# Top
	draw color $color1
        draw triangle $outer_1 $outer_2 $x2
        draw triangle $outer_1 $outer_8 $x2
        draw triangle $outer_5 $outer_4 $x2
        draw triangle $outer_5 $outer_6 $x2
	draw color $color2
        draw triangle $outer_3 $outer_2 $x2
        draw triangle $outer_3 $outer_4 $x2
        draw triangle $outer_7 $outer_6 $x2
        draw triangle $outer_7 $outer_8 $x2

}
