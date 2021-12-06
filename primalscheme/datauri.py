# This module exports functions to convert text or binary data to a data URI readable by
# igv.js.
# ALL CREDIT: https://github.com/igvteam/igv-reports by https://github.com/jrobinso

from typing import Union

from base64 import b64encode
from gzip import compress


def get_data_uri(data: Union[str, bytes]) -> str:

    """
    Return a data uri for the input, which can be either a string or byte array
    """

    if isinstance(data, str):
        data = compress(data.encode())
        mediatype = "data:application/gzip"
    else:
        if data[0] == 0x1F and data[1] == 0x8B:
            mediatype = "data:application/gzip"
        else:
            mediatype = "data:application:octet-stream"

    enc_str = b64encode(data)

    data_uri = mediatype + ";base64," + str(enc_str)[2:-1]
    return data_uri
