#! vmd -dispdev -eof < g2.tcl >& log

#change all the pdb to the same coordinate system here we set index 0 point as the origin
#and index 3 in y axis , index 1 in xy plane
#then we could calculate SDF, we should calculate the box size and grids
#it will print sdf in cube format
#rewrite cube format by yourself, such as add molecular number, box origin point, grid number and coordinate (unit in bohr)  
#
#=====================================
# spatial distribution function
# @ hzau.wuhan.cn/
# @ liufang
# @
#2016-10-20
#=====================================

#sleep N second
proc sleep {N} {
    after [expr {int($N*1000)}]
}

proc transcoord {file} {
    #https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    #http://www.cnblogs.com/graphics/archive/2012/08/08/2609005.html
    #http://blog.sina.com.cn/s/blog_4a4a8c7d0100o4hq.html
    puts $file 
    mol new $file type pdb

    set all [atomselect top all]
    set aatom [atomselect top "index 0"]
    set gc [veczero]

    set origin [lindex [$aatom get {x y z}] 0]
    set origin [vecsub $gc $origin]
    $all moveby $origin
    #$all writepdb 11.pdb
    #now index 0 atom coordinate is 0.0 0.0 0.0
    #then translate theta with y axis
    #translate around yaxis
    set batom [atomselect top "index 3"]
    set bc [lindex [$batom get {x y z}] 0]
    set bcy [lindex $bc 2]
    puts $bc
    set bc [lreplace $bc 1 1 0.0]
    set xvec {1.0 0.0 0.0}
    set xb [vecdot $bc $xvec]
    puts $xb
    set bclength [veclength $bc]
    set cost [expr $xb/$bclength]
    puts $cost
    ####
    #rad to angle
    set PI [expr 2*asin(1.0)]
    ####

    set xtheta [expr (acos($cost))*180/$PI]

    #set matrix [transaxis y -152.382]
    #$all move $matrix
    #$all writepdb 12.pdb
    if { $bcy < 0.0 } {
        set ymatrix [transaxis y [expr 0-$xtheta]]
        $all move $ymatrix
    } else {
        set ymatrix [transaxis y $xtheta]
        $all move $ymatrix
    }
    #$all writepdb 113.pdb

    #translate around zaxis
    set batom [atomselect top "index 3"]
    set bc [lindex [$batom get {x y z}] 0]
    set bcz [lindex $bc 0]
    set yvec {0.0 1.0 0.0}
    set yb [vecdot $bc $yvec]
    puts $yb
    set bclength [veclength $bc]
    set cosa [expr $yb/$bclength]
    puts $cosa
    ####
    #rad to angle
    ####
    set yalpha [expr (acos($cosa))*180/$PI]
    if { $bcz < 0.0 } {
        set zmatrix [transaxis z [expr 0-$yalpha]]
        $all move $zmatrix
    } else {
        set zmatrix [transaxis z $yalpha]
        $all move $zmatrix
    }
    #$all writepdb 114.pdb

    #translate around yaxis
    set catom [atomselect top "index 1"]
    set cc [lindex [$catom get {x y z}] 0]
    puts $cc
    set ccz [lindex $cc 2]
    set cc [lreplace $cc 1 1 0.0]
    set xvec {1.0 0.0 0.0}
    set cx [vecdot $cc $xvec]
    puts $cx
    set cclength [veclength $cc]
    set cosb [expr $cx/$cclength]
    ####
    #rad to angle
    ####
    set xbeta [expr (acos($cosb))*180/$PI]
    puts $xbeta
    if { $ccz < 0.0 } {
        set ymatrix [transaxis y [expr 0-$xbeta]]
        $all move $ymatrix
    } else {
        set ymatrix [transaxis y $xbeta]
        $all move $ymatrix
    }
    $all writepdb $file.pdb
    
mol delete top
}

proc tgrep { pattern filename} {
    global data
    set file [open $filename r]
    foreach line [split [read $file] \n] {
	if { [regexp $pattern $line]} {
	#puts $line
        set str [string trim $line]
	lappend data $str
	#puts "haha "
	#}
    }
    close $file
return $data
}
####
#main function
set path "."
set files [glob "$path/*"]
set numpdb 0
foreach f $files {
    if {[string match *.pdb $f]} {
        set numpdb [expr $numpdb + 1]
        #puts $f
        set file [lindex [split $f /] 1]
        #puts $file
        lappend newpdb $file.pdb
	transcoord $file
        }
    sleep 0.01
}
#puts $newpdb
#puts $numpdb
sleep 10
####
#read newpdb
foreach pdb $newpdb {
    set pattern {O1  MOH}
    #this proc just like grep command of UNIX
    puts $pdb
    tgrep $pattern $pdb
    ####
    #heihei now all the select atom which will be analyzed is collected in data
    #puts $data
    



}
puts $data
####
#all the space were set to 0
#origin coord come from vmd measure center, maybe it is usefull!
#set originx -2
#set originy 3
#set originz 1
#minior coord
#0.1A is the length of each step,and i have measured the box with vmd,here each step is 0.1 A
set step 0.1
set numx 150
set numy 200
set numz 150
#now the minior point is -9.5 -7 -6.5
set minx -9.5
set miny -7.0
set minz -6.5
for {set i 0} {$i<=$numx} {incr i 1} {
    for {set j 0} {$j<=$numy} {incr j 1} {
        for {set k 0} {$k<=$numz} {incr k 1} {
	    set coord 0
            set array($i,$j,$k) $coord
       	    #puts $array(-9.5,-6.0,-6.0)
	    #puts array($coordx,$coordy,$coordz)
            #set array(-0.0,-0.0,-0.0) 0
        }
    }
}

foreach line $data {
    set lline [llength $line]
    if {$lline == 11} {
        set pointx [lindex $line 6]
	set pointy [lindex $line 7]
	set pointz [lindex $line 8]
    } else {
        set pointx [lindex $line 5]
	set pointy [lindex $line 6]
	set pointz [lindex $line 7]
    }
    set x [expr int(($pointx-$minx)/$step)]
    set y [expr int(($pointy-$miny)/$step)]
    set z [expr int(($pointz-$minz)/$step)]
    #set pointx [format "%.1f" [expr $minx+$x*$step]]
    #set pointy [format "%.1f" [expr $miny+$y*$step]]
    #set pointz [format "%.1f" [expr $minz+$z*$step]]
    #if {$pointx == -0.0} {
    #set pointx 0.0
    #}
    #if {$pointy == -0.0} {
    #set pointy 0.0
    #}
    #if {$pointz == -0.0} {
    #set pointz 0.0
    #}
    set array($x,$y,$z) [incr array($x,$y,$z) 1]
    puts array($x,$y,$z)
    puts $array($x,$y,$z)

}
set numdata [llength $data]
set cubefile sdf.cube
set file [open $cubefile w]
for {set i 0} {$i<=$numx} {incr i 1} {
    for {set j 0} {$j<=$numy} {incr j 1} {
        for {set k 0} {$k<=$numz} {incr k 1} {
	    #set pointx [format "%.1f" [expr ($minx+$i*$step)/$numdata]]
	    #set pointy [format "%.1f" [expr ($miny+$j*$step)/$numdata]]
	    #set pointz [format "%.1f" [expr ($minz+$k*$step)/$numdata]]
	    #if {$pointx == -0.0} {
	    #set pointx 0.0
	    #}
	    #if {$pointy == -0.0} {
	    #set pointy 0.0
	    #}
	    #if {$pointz == -0.0} {
	    #set pointz 0.0
            #}    
	    puts -nonewline $file "$array($i,$j,$k)\ "
	    if {[expr $k%6]==5} {
	        puts -nonewline $file "\n"
	    }
	    }
	#puts $array()
	puts -nonewline $file "\n"
    }
}
close $file
