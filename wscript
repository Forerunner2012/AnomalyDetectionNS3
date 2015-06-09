# -*- Mode: python; py-indent-offset: 4; indent-tabs-mode: nil; coding: utf-8; -*-

# def options(opt):
#     pass

# def configure(conf):
#     conf.check_nonfatal(header_name='stdint.h', define_name='HAVE_STDINT_H')

def build(bld):
    module = bld.create_ns3_module('ad-detector-simulator', ['core'])
    module.source = [
        'model/ad-detector-simulator.cc',
        'helper/ad-detector-simulator-helper.cc',
        ]

    module_test = bld.create_ns3_module_test_library('ad-detector-simulator')
    module_test.source = [
        'test/ad-detector-simulator-test-suite.cc',
        ]

    headers = bld(features='ns3header')
    headers.module = 'ad-detector-simulator'
    headers.source = [
        'model/ad-detector-simulator.h',
        'helper/ad-detector-simulator-helper.h',
        ]

    if bld.env.ENABLE_EXAMPLES:
        bld.recurse('examples')

    # bld.ns3_python_bindings()

