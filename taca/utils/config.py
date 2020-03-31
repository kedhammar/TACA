""" Load and parse configuration file
"""
import yaml

CONFIG = {}

def load_config(config_file):
    """Loads a configuration file.
    """
    config = {}
    if type(config_file) is file:
        config.update(yaml.load(config_file, Loader=yaml.FullLoader) or {})
        return config
    else:
        try:
            with open(config_file, 'r') as f:
                content = yaml.load(f, Loader=yaml.FullLoader)
                config.update(content)
                return content
        except IOError as e:
            e.message = "Could not open configuration file \"{}\".".format(config_file)
            raise e



def load_yaml_config(config_file):
    """Load YAML config file

    :param str config_file: The path to the configuration file.

    :returns: A dict of the parsed config file.
    :rtype: dict
    :raises IOError: If the config file cannot be opened.
    """
    if type(config_file) is file:
        CONFIG.update(yaml.load(config_file, Loader=yaml.FullLoader) or {})
        return CONFIG
    else:
        try:
            with open(config_file, 'r') as f:
                content = yaml.load(f, Loader=yaml.FullLoader)
                CONFIG.update(content)
                return content
        except IOError as e:
            e.message = "Could not open configuration file \"{}\".".format(config_file)
            raise e
