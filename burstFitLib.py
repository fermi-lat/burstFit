def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['burstFit'], package = 'burstFit')
    env.Tool('evtbinLib')
    env.Tool('optimizersLib')
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_graphLib')
    env.Tool('st_streamLib')
    env.Tool('tipLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    

def exists(env):
    return 1
