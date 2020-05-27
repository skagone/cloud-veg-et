import argparse

from vegetLib.veget import VegET

from vegetLib.log_logger import log_make_logger

def get_parser():
    parser = argparse.ArgumentParser(description='Run the veget model with specified tile')
    parser.add_argument('tile', metavar='TILE', type=str, nargs='*',
            help='the tile to process - example: 40N-80E')
    parser.add_argument('-c', '--configdir', help='specify and alternate config dir example: -c sample_config ', default='./sample_config', type=str)
    return parser


def command_line_runner():
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['configdir']:
        print("configdir", args['configdir'])

    if args['tile']:
        print("tile", args['tile'])

    tile = args['tile'][0]

    config_directory = args['configdir']

    print('Using config dir {} to process tile {}'.format(config_directory, tile))

    log.info('Processing Tile {}'.format(tile))

    log.error('Someday tony needs to refine logging and determine a std format for ELK')
    log.error('or logging agents and logging backends ... docker deployments')

    # RUN the class Veget
    myveg = VegET(config_directory)

if __name__ == '__main__':
    log = log_make_logger('CLOUD_VEGET')
    command_line_runner()
