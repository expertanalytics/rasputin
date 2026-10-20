"""Signed download links: HMAC-SHA256 over ``user/name/filename/expires``, in hex.

Pure: no I/O and no clock; ``now`` is passed in. A link opens one stored file
until ``expires`` (unix seconds) and is checked in constant time.
"""

from __future__ import annotations

import hashlib
import hmac

__all__ = ["sign", "url", "valid"]


def sign(key: bytes, user: str, name: str, filename: str, expires: int) -> str:
    """The signature of one link, as lowercase hex."""
    message = f"{user}/{name}/{filename}/{expires}".encode()
    return hmac.new(key, message, hashlib.sha256).hexdigest()


def valid(
    key: bytes, user: str, name: str, filename: str, expires: int, sig: str, now: float
) -> bool:
    """Whether ``sig`` signs these parts and ``now`` is before ``expires``."""
    good = hmac.compare_digest(sign(key, user, name, filename, expires).encode(), sig.encode())
    return good and now < expires


def url(public_url: str, user: str, name: str, filename: str, expires: int, sig: str) -> str:
    """The link on ``public_url``: ``/files/<user>/<name>?fn=&exp=&sig=``; every part is
    of URL-safe characters (the route refuses any other), so none is escaped."""
    query = f"fn={filename}&exp={expires}&sig={sig}"
    return f"{public_url.rstrip('/')}/files/{user}/{name}?{query}"
