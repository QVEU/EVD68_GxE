#!/usr/bin/env python

"""
Adam Catching, Ph. D.
2022-08-29
NIH-NIAID-LVD
Quantitative Virology and Evolution Unit
"""

import sys
import os
from opentrons import protocol_api

metadata = {
    'protocolName': 'Viral Passaging, New Passage',
    'author': 'Adam Catching',
    'description': 'Viral Passaging, Serially Seed',
    'apiLevel': '2.9'
}

def run(ctx: protocol_api.ProtocolContext):

    """
    ============================================================================
                                     Set-up
    ============================================================================
    """

    # Load tips
    tiprack1 = ctx.load_labware("opentrons_96_filtertiprack_20ul", '11')
    tiprack2 = ctx.load_labware("opentrons_96_filtertiprack_20ul", '9')
    # Load plate locations
    plate1 = ctx.load_labware("corning_96_wellplate_360ul_flat", '7')
    plate2 = ctx.load_labware("corning_96_wellplate_360ul_flat", '4')
    plate3 = ctx.load_labware("corning_96_wellplate_360ul_flat", '8')
    plate4 = ctx.load_labware("corning_96_wellplate_360ul_flat", '5')
    plates = [plate1, plate2, plate3, plate4]

    # Plate select samples are transfered to
    tempdeck = ctx.load_module('Temperature Module Gen2', '10')
    tempdeck.set_temperature(8)
    final_plate = tempdeck.load_labware(
                'corning_96_wellplate_360ul_flat',
                'consolidation plate')

    # Load multi-channel pipette
    m20 = ctx.load_instrument('p20_multi_gen2', mount='left',
    tip_racks=[tiprack1, tiprack2])

    """
    ============================================================================
                                      Run
    ============================================================================
    """

    # Define row prefixes
    rows = [str(x) for x in 'ABCDEFGH']
    # Iterate over the plate locations, i is the plate/final column number
    for i in range(3):
        for j, plate in enumerate(plates):
            m20.transfer(10,
            final_plate.wells_by_name()[''.join(['A', str(j+1+(4 * i))])],
            plate.wells_by_name()['A1'], replace_tips='never'
            )
            # Iterate over rows, j is the final row for the samples

            # Define old location from .csv and new location by which plate
            old_loc = [''.join(['A', str(x+1)]) for x in range(11)]
            new_loc = [''.join(['A', str(x+2)]) for x in range(11)]
            # Transfer from specific condition plate to consolidated plate
            m20.transfer(10,
            plate.wells('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11'),
            plate.wells('A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12'),
            mix_before=(3, 5), replace_tips='never')

            ctx.delay(minutes=1, msg=f'Replace plates')
