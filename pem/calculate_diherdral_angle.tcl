#! /usr/bin/env tcl

# vmd embeded

proc cal_dihedral {file} {
    mol new $file type xyz
    global i
    global flog
    global a
    global b
    global c
    global d
    global dih
    set D [measure dihed "$a $b $c $d"]
    set delta_dih [expr $D-$dih]
    puts -nonewline $flog [format "%6d" $i]
    puts -nonewline $flog [format "%12.5f" $D]
    puts -nonewline $flog [format "%12.5f" $delta_dih]
    puts -nonewline $flog "\n"
    mol delete top
}

#main function
#read min.dat file
set dat "min.dat"
set data ""
set inStream [open $dat r]

foreach line [split [read $inStream] \n] {
    set line [string trim $line]
    set rec [split $line " "]
    set num [llength $rec]
    if {$num<=5} {
        set traj [lindex $line 0]
        puts $traj
    } else {
        for {set i 0} {$i<$num} {incr i 1} {
            set geom [lindex $rec $i]
            puts $geom
            lappend data $traj-$geom
        }
    }
    
}
#puts $data
#read geom file
set cal_file cal.dat
set flog [open $cal_file w]
set i 0

set a 3
set b 12
set c 13
set d 15
#calculate initial geom
set initial_file geom.xyz
mol new $initial_file type xyz
####diheral
set dih [measure dihed "$a $b $c $d"]
foreach geom $data {
    set path "../geom/"
    puts $geom
    set geom_file $path$geom
    cal_dihedral $geom_file
    incr i 1
}
close $flog
