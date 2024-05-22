#!/usr/bin/python

"""
Adam Catching, Ph. D.
2023-03-21
NIH-NIAID-LVD
Quantitative Virology and Evolution Unit
"""


from opentrons import protocol_api

metadata = {
    'protocolName': 'TCID50, 4 technical triplicates',
    'author': 'Adam Catching',
    'description': 'TCID50, in triplicate, of 4 samples',
    'apiLevel': '2.9'
}

def run(ctx: protocol_api.ProtocolContext):

    """
    ============================================================================
                                     Set-up
    ============================================================================
    """


    # Load tips
    tiprack1 = ctx.load_labware("opentrons_96_filtertiprack_200ul", '7')
    tiprack2 = ctx.load_labware("opentrons_96_filtertiprack_200ul", '8')

    # Load plate locations
    plate1 = ctx.load_labware("corning_96_wellplate_360ul_flat", '4')
    plate2 = ctx.load_labware("corning_96_wellplate_360ul_flat", '5')
    plate3 = ctx.load_labware("corning_96_wellplate_360ul_flat", '6')
    plate4 = ctx.load_labware("corning_96_wellplate_360ul_flat", '1')
    plate5 = ctx.load_labware("corning_96_wellplate_360ul_flat", '2')
    plate6 = ctx.load_labware("corning_96_wellplate_360ul_flat", '3')

    # Define serial dilution plate
    wellplate = ctx.load_labware('usascientific_96_wellplate_2.4ml_deep', '11')

    # Load sample locations
    tubes = 'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap'
    tuberack = ctx.load_labware(tubes, '10')

    # Load multi-channel pipette
    p300 = ctx.load_instrument('p300_single_gen2', mount='right',
    tip_racks=[tiprack1])
    # Load single channel pipette
    m300 = ctx.load_instrument('p300_multi_gen2', mount='left',
    tip_racks=[tiprack2])

    # Name samples
    sample1 = tuberack.wells_by_name()["A1"]
    sample2 = tuberack.wells_by_name()["A2"]
    sample3 = tuberack.wells_by_name()["A3"]
    sample4 = tuberack.wells_by_name()["A4"]

    samples = [sample1, sample2, sample3, sample4]

    # Define the initial locations of the serial dilution as a list
    init_locs = ["".join(['A', str(i+1)]) for i in range(12)]
    initial_dil = [wellplate.wells_by_name()[x] for x in init_locs]


    """
    ============================================================================
                                     Protocol
    ============================================================================
    """

    # Pick up tip and move 166 uL from sample A1 to multiwell A1
    def transfer_sample(start_loc, end_loc):
        p300.pick_up_tip()
        p300.mix(4, 100, start_loc)
        p300.aspirate(166, start_loc)
        p300.dispense(166, end_loc)
        p300.drop_tip()

    # Transfer samples A1 through A4 to deepwell samples
    for i, sample in enumerate(samples):
        for j in range(3):
            # Transfer stock sample into 2 consecutive deepwells
            p300.transfer(166, sample.bottom(5),
            initial_dil[i * 3 + j].bottom(10), blow_out=True,
            air_gap=True, mix_before=(2, 100))
            #transfer_sample(sample, initial_dil[i * 3 + j])

    # Serially dilute row by row
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for i, old_row in enumerate(rows[:-1]):
        # Define new row of undiluted media
        new_row = rows[i+1]
        old_locs = ["".join([old_row, str(i+1)]) for i in range(12)]
        new_locs = ["".join([new_row, str(i+1)]) for i in range(12)]

        # Define the rows that have (old_source) and have not (new_source) been
        # serially diluted into.
        old_source = [wellplate.wells_by_name()[x].bottom(10) for x in old_locs]
        new_source = [wellplate.wells_by_name()[x].bottom(10) for x in new_locs]

        # Transfer to next 10-fold dilution (A1 -> B1)
        p300.transfer(166.6, old_source,new_source,
        mix_before=(3, 100),new_tip='always')

    # Transfer serial diluted virus to the seeded plates of cells
    m300.transfer(111.1, wellplate.wells_by_name()['A1'].bottom(3),
    plate1.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A2'].bottom(3),
    plate2.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A3'].bottom(3),
    plate3.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A4'].bottom(3),
    plate4.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A5'].bottom(3),
    plate5.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A6'].bottom(3),
    plate6.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

     # All displayed 96-well plates are titered, replace with the next batch of 6
    ctx.delay(minutes=1, msg=f'Replace plates')

    m300.transfer(111.1, wellplate.wells_by_name()['A7'].bottom(3),
    plate1.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A8'].bottom(3),
    plate2.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)
    
    m300.transfer(111.1, wellplate.wells_by_name()['A9'].bottom(3),
    plate3.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)
    
    m300.transfer(111.1, wellplate.wells_by_name()['A10'].bottom(3),
    plate4.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A11'].bottom(3),
    plate5.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

    m300.transfer(111.1, wellplate.wells_by_name()['A12'].bottom(3),
    plate6.columns('1', '2', '3','4', '5', '6', '7',
    '8', '9', '10', '11', '12'),
    mix_before=(3, 100), air_gap=20)

