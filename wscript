#! python


def build(ctx):

 ctx(
        features='run_m_script',
        source='estimation.m',
        deps=ctx.path_to(ctx, 'OUT_DATA', 'bullionist_controversy.mat'),
        target=[
			ctx.path_to(ctx, 'OUT_FIGURES', 'impulse_response.pdf'),
			
			ctx.path_to(ctx, 'OUT_FIGURES', 'fevd.pdf' )],
        name='estimation'
    )
