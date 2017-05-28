#! vmd -dispdev -eof < *.tcl >& log
#read dcd file and calculate hbond number of selection1 and selection2 and select these resid 
# then write pdb file in different file as different selection2 residue number
#the last line of log file will print hbond resdiue number of selection2 which interact with selection1
# and from these selection2 number we can understand the distribution of hbond number between the two selection

#=====================================
# spatial distribution function
# @ hzau.wuhan.cn/
# @ liufang
# @
#2016-10-20
#=====================================


mol new {./complex.pdb} type pdb
mol off top
set output_frame_offset 0
set first_frame 0
set last_frame -1
set filelist {./sin_equ.dcd}
# set outputPath {};	# Output path for results files. If it is empty, output is written to the same location as an input file

#Set selection to look for atoms

set selection1 [atomselect top "residue 0 and name H H1 H3 H4 H5 H6 H7 H8 H9 H10 H11 O1 O2 O3 O4"]
set selection2 [atomselect top "resname MOH and name HC1 HC2 HC3 HO1 O1"]


#Set H-bond definition criteria
set cutoff 3.81
set angle 45.1

#load file
animate read dcd $filelist beg $first_frame end $last_frame waitfor all

#Extract frames from file
set num_steps [molinfo top get numframes]
#Open files for writing
set filename "1.dat"	
#set fid [open $filename w]
#set fidlog [open $filename.log w]

#puts $fidlog "Coordinates from file $filelist read"
#puts $fidlog "Output file $filename initiated"
	
set acceptor "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
set donor "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
set all "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
set acceptor_n [split $acceptor]
set donor_n [split $donor]
set all_n [split $all]

set res_file ""
for {set frame 0} {$frame < $num_steps} {incr frame} {
    #Update the frames
    $selection1 frame $frame
    $selection1 update

    $selection2 frame $frame
    $selection2 update

    #Count H-bonds in selection
    set ceramide_water_hBonds [measure hBonds $cutoff $angle $selection1 $selection2]
    set water_ceramide_hBonds [measure hBonds $cutoff $angle $selection2 $selection1]
#    puts $fidlog "H-bonds measured for frame $frame"
    #puts $ceramide_water_hBonds	
    #puts "123"
    #puts $water_ceramide_hBonds
    #puts "456"
    set ceramide_water_hBonds_full {}
    puts "123"
    puts $frame
    #puts $water_ceramide_hBonds
    set res ""
    ###
    #residue number hbond with residue 0
    set res_s ""
    ###
    foreach donor [lindex $ceramide_water_hBonds 0] acceptor [lindex $ceramide_water_hBonds 1] {
        #puts $donor; #puts $acceptor
        set sel [atomselect top "index $donor"]
        set donorResidue [$sel get residue]
	#puts $donorResidue
        lappend res $donorResidue
        
        $sel delete
			
	set sel [atomselect top "index $acceptor"]
	set acceptorResidue [$sel get residue]
        #puts $sel; #puts $acceptorResidue
        lappend res $acceptorResidue

       	if {$acceptorResidue ni $res_s} {
            lappend res_s $acceptorResidue
        }

        $sel delete
			
	lappend ceramide_water_hBonds_full "$acceptorResidue $donorResidue a"
    }
    set ceramide_water_hBonds_full [lsort -integer -index 0 $ceramide_water_hBonds_full]
#    puts $fidlog "_________________ ceramide_water_hBonds_full $ceramide_water_hBonds_full"
    puts $ceramide_water_hBonds_full		
    set water_ceramide_hBonds_full {}
    foreach donor [lindex $water_ceramide_hBonds 0] acceptor [lindex $water_ceramide_hBonds 1] {
	set sel [atomselect top "index $donor"]
	set donorResidue [$sel get residue]
	lappend res $donorResidue
        if {$donorResidue ni $res_s} {
            lappend res_s $donorResidue
        }

        $sel delete
			
	set sel [atomselect top "index $acceptor"]
	set acceptorResidue [$sel get residue]
	lappend res $acceptorResidue
        $sel delete
			
	lappend water_ceramide_hBonds_full "$donorResidue $acceptorResidue d"
    }

    set water_ceramide_hBonds_full [lsort -integer -index 0 $water_ceramide_hBonds_full]
    #puts $fidlog "___"
    set all_hBonds [lsort -integer -index 0 [concat $ceramide_water_hBonds_full $water_ceramide_hBonds_full]];   # Array with All the hydrogen bonds water makes with ceramde
    puts $water_ceramide_hBonds_full
# 	puts $fidlog "An array with all the H-bonds created for frame $frame"
    set all_hBonds_number [llength $all_hBonds]
    #puts $fidlog "all_hBonds_number for frame $frame: $all_hBonds_number"
    
    #count residue number for each frame
    set A [llength $res_s]   
    set A_new [expr [lindex $all_n $A]+1]
    set all_n [lreplace $all_n $A $A $A_new]
    puts $all_n

    puts $res
    puts $res_s

    ####
    #sin and interact moh
    if {$A in $res_file} {
        cd $A
        set atom [atomselect top "residue $res_s 0" frame $frame]
        $atom writepdb $frame.pdb
        cd ..
    } else {
        file mkdir $A
        cd $A
        set atom [atomselect top "residue $res_s 0" frame $frame]
        $atom writepdb $frame.pdb
        cd ..
        ###
    }

    #count hbond number for each frame
    #set a [llength $ceramide_water_hBonds_full]
    #set a_new [expr [lindex $acceptor_n $a]+1]
    
    #set acceptor_n [lreplace $acceptor_n $a $a $a_new]
    #
    #set d [llength $water_ceramide_hBonds_full]
    #set d_new [expr [lindex $donor_n $d]+1]
    #set donor_n [lreplace $donor_n $d $d $d_new]
    #
    #set A [llength $all_hBonds]
    #set A_new [expr [lindex $all_n $A]+1]
    #set all_n [lreplace $all_n $A $A $A_new]
    #puts $adnumberlog "acceptor and donor and all_hBonds_number for frame $frame: "
    
    #puts hbond and residue 0 coord
    #puts all the pdb 
    #set all [atomselect top all frame $frame]
    #$all writepdb $frame.pdb
}
###
#count_hbond residue number for all frame
puts $all_n
#count hbond number for all frame
#puts $acceptor_n
#puts $donor_n
#puts $all_n
###

mol on top


#bell
puts "Finished !!!"




