## -*- Mode: python; py-indent-offset: 4; indent-tabs-mode: nil; coding: utf-8; -*-

def build(bld):
    obj = bld.create_ns3_program('ads-simulator',
                                 ['core'])
    obj.source = 'ads-simulator.cc'

    obj = bld.create_ns3_program('wave-simulator',
                                 ['core', 'applications', 'mobility', 'network', 'wifi','wave'])
    obj.source = 'wave-simulator.cc'
    
    obj = bld.create_ns3_program('tranform_rawdata_to_ns3mobilityhelper.cc',
                                 ['core', 'applications', 'mobility', 'network', 'wifi','wave'])
    obj.source = 'tranform_rawdata_to_ns3mobilityhelper.cc'
    


