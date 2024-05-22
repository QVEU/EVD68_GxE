#!/usr/bin/env python

"""
Adam Catching, Ph. D.
2022-08-26
NIH-NIAID-LVD
Quantitative Virology and Evolution Unit
"""

from opentrons import protocol_api

metadata = {
    'protocolName': 'Sample Transfer',
    'author': 'Adam Catching',
    'description': 'Sample Transfer, Passaing',
    'apiLevel': '2.9'
}

def run(ctx: protocol_api.ProtocolContext):

    """
    ============================================================================
                                     Set-up
    ============================================================================
    """

    # Load tips
    tiprack1 = ctx.load_labware("opentrons_96_filtertiprack_200ul", '11')
    # Load plate locations
    plate1 = ctx.load_labware("corning_96_wellplate_360ul_flat", '7')
    plate2 = ctx.load_labware("corning_96_wellplate_360ul_flat", '4')
    plate3 = ctx.load_labware("corning_96_wellplate_360ul_flat", '1')
    plate4 = ctx.load_labware("corning_96_wellplate_360ul_flat", '9')
    plate5 = ctx.load_labware("corning_96_wellplate_360ul_flat", '6')
    plate6 = ctx.load_labware("corning_96_wellplate_360ul_flat", '3')
    plates = [plate1, plate2, plate3, plate4, plate5, plate6]
    # Plate select samples are transfered to
    tempdeck = ctx.load_module('Temperature Module Gen2', '10')
    tempdeck.set_temperature(8)
    final_plate = tempdeck.load_labware('corning_96_wellplate_360ul_flat','elution plate')
    # Load multi-channel pipette
    p300 = ctx.load_instrument('p300_single_gen2', mount='right', tip_racks=[tiprack1])

    """
    ============================================================================
                                      Run
    ============================================================================
    """
    # Array of locations needs to be read in then transposed
    # No current way for OpenTrons on Windows to read a .csv file, has to be
    # manually entered
    loc_data = ['0,6,0,0,0,0,9,12,0,12,8,12', '2,8,0,12,0,0,0,12,0,0,0,12',
    '2,0,0,12,0,8,0,12,0,11,0,12', '3,0,0,12,0,0,6,12,0,12,0,12',
    '3,4,0,12,0,12,0,12,0,12,1,12', '2,8,0,12,0,0,0,12,0,0,0,12',
    '2,12,0,12,0,7,0,12,0,12,0,12', '4,12,0,12,0,12,0,12,0,0,0,12']
    # Initiate location array
    locations = []
    #f.close()
    # Iterate over cleaned up lines of .csv file
    for i, locs in enumerate(loc_data):
        # Define new row
        location_row = []
        # Split csv row into list to be iterated over
        for loc in locs.split(','):
            # Convert possible string value to integer
            int_loc = int(loc)
            # Make sure the plate locations are now starting at 1 and
            # values too low are set at a minimum of 1
            if int_loc < 5:
                MOI_1_loc = 1
            else:
                MOI_1_loc = int_loc - 4
            # Append row of plate locations to array
            location_row.append(MOI_1_loc)
        locations.append(location_row)
    print(locations)
    # Define transposing function
    T = lambda A: [[A[l][m] for l in range(len(A))] for m in range(len(A[0]))]
    # Transpose locations to match 1 row to one plate
    locations_T = T(locations)
    #print(locations_T)
    # Define row prefixes
    rows = [str(x) for x in 'ABCDEFGH']
    # Iterate over the plate locations, i is the plate/final column number
    for i, plate in enumerate(plates):
        # Iterate over rows, j is the final row for the samples
        for j, row in enumerate(rows):
            # Define old location from .csv and new location by which plate
            old_loc = ''.join([row, str(locations_T[i][j])])
            new_loc = ''.join([row, str(i+1)])
            # Transfer from specific condition plate to consolidated plate
            p300.transfer(100, plate.wells_by_name()[old_loc].bottom(1),
            final_plate.wells_by_name()[new_loc].bottom(2),
            mix_before=(3, 50),new_tip='always')

    ctx.delay(minutes=1, msg=f'Replace plates\n\n\n')

    # Iterate over the plate locations, i is the plate/final column number
    for i, plate in enumerate(plates):
        # Iterate over rows, j is the final row for the samples
        for j, row in enumerate(rows):
            # Define old location from .csv and new location by which plate
            old_loc = ''.join([row, str(locations_T[i+6][j])])
            new_loc = ''.join([row, str(i+7)])
            # Transfer from specific condition plate to consolidated plate
            p300.transfer(100, plate.wells_by_name()[old_loc].bottom(1),
            final_plate.wells_by_name()[new_loc].bottom(2),
            mix_before=(3, 50),new_tip='always')
