#!/usr/bin/env python3
"""
transit_finder.py - find ISS and Tiangong transits of the Sun and Moon.

Scans the next N days for passes where a space station crosses the solar or
lunar disc as seen from one observing site. Reports true disc crossings, plus
"possible" near misses far enough out that orbit-prediction drift could still
move the path over the site.

The observing site is an input only. It is never written to the repo:
  1. --lat / --lon / --elev on the command line, or
  2. TT_LAT / TT_LON / TT_ELEV environment variables, or
  3. a local location.json (gitignored), see location.example.json.

Orbit data (TLEs), first one found wins:
  1. --tle FILE
  2. data/stations.tle in this repo (refreshed daily by the GitHub Action)
  3. a live download from CelesTrak
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import urllib.request
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from skyfield.api import EarthSatellite, Loader, load, wgs84

HERE = Path(__file__).resolve().parent

STATIONS = {25544: "ISS", 48274: "Tiangong"}
CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR={catnr}&FORMAT=TLE"

SUN_RADIUS_KM = 696_000.0
MOON_RADIUS_KM = 1_737.4

COARSE_STEP_S = 0.5     # sampling step across a pass
FINE_STEP_S = 0.002     # refinement step around closest approach


# ---------------------------------------------------------------- inputs

def load_location(args) -> tuple[float, float, float]:
    lat, lon, elev = args.lat, args.lon, args.elev
    if lat is None or lon is None:
        if os.environ.get("TT_LAT") and os.environ.get("TT_LON"):
            lat = float(os.environ["TT_LAT"])
            lon = float(os.environ["TT_LON"])
            elev = float(os.environ.get("TT_ELEV", elev or 0.0))
        else:
            cfg = Path(args.location_file)
            if cfg.is_file():
                data = json.loads(cfg.read_text())
                lat, lon = float(data["lat"]), float(data["lon"])
                elev = float(data.get("elev_m", elev or 0.0))
    if lat is None or lon is None:
        sys.exit("No observing site given. Use --lat/--lon, TT_LAT/TT_LON, "
                 "or a location.json (see location.example.json).")
    return lat, lon, float(elev or 0.0)


def parse_tle_text(text: str) -> dict[int, tuple[str, str, str]]:
    """Return {catnr: (name, line1, line2)} for every 3-line or 2-line TLE."""
    lines = [ln.rstrip() for ln in text.replace("\r", "").splitlines() if ln.strip()]
    out: dict[int, tuple[str, str, str]] = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name, l1, l2 = "", lines[i], lines[i + 1]
            i += 2
        elif (i + 2 < len(lines) and lines[i + 1].startswith("1 ")
              and lines[i + 2].startswith("2 ")):
            name, l1, l2 = lines[i].strip(), lines[i + 1], lines[i + 2]
            i += 3
        else:
            i += 1
            continue
        out[int(l1[2:7])] = (name, l1, l2)
    return out


def load_tles(args) -> tuple[dict[int, tuple[str, str, str]], str]:
    candidates = [Path(args.tle)] if args.tle else [HERE / "data" / "stations.tle"]
    for path in candidates:
        if path.is_file():
            tles = parse_tle_text(path.read_text())
            if all(c in tles for c in STATIONS):
                return tles, str(path)
            if args.tle:
                sys.exit(f"{path} is missing one of {sorted(STATIONS)}")
    text = ""
    try:
        for catnr in STATIONS:
            req = urllib.request.Request(CELESTRAK_URL.format(catnr=catnr),
                                         headers={"User-Agent": "transit-finder"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                text += resp.read().decode() + "\n"
    except OSError as exc:
        sys.exit(f"No local TLE file and CelesTrak download failed: {exc}")
    tles = parse_tle_text(text)
    if not all(c in tles for c in STATIONS):
        sys.exit("CelesTrak returned no usable TLE for one of the stations.")
    return tles, "CelesTrak (live)"


def load_ephemeris():
    """de421 from the skyfield-data package if installed, else download once."""
    try:
        import warnings
        with warnings.catch_warnings():
            # skyfield-data warns about its bundled Earth-orientation file;
            # only its de421.bsp is used here, so the warning doesn't apply.
            warnings.simplefilter("ignore")
            import skyfield_data  # type: ignore
            path = Path(skyfield_data.get_skyfield_data_path()) / "de421.bsp"
        if path.is_file():
            return Loader(str(path.parent), verbose=False)("de421.bsp")
    except ImportError:
        pass
    cache = HERE / "cache"
    cache.mkdir(exist_ok=True)
    return Loader(str(cache), verbose=False)("de421.bsp")


# ---------------------------------------------------------------- geometry

@dataclass
class Event:
    station: str
    body: str
    kind: str                 # "transit" or "possible"
    time_utc: str
    time_local: str
    days_out: float
    min_sep_arcmin: float
    disc_radius_arcmin: float
    duration_s: float         # time inside the disc (0 for a near miss)
    station_alt_deg: float
    station_az_deg: float
    station_range_km: float
    body_alt_deg: float
    body_az_deg: float
    direction: str            # compass direction of the event, e.g. "SE"
    low: bool                 # under --low-alt: may be hidden by trees/houses
    offset_km: float          # how far the site sits from the path centerline
    path_halfwidth_km: float  # half-width of the ground strip that sees a transit
    tle_age_days: float


COMPASS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]


def compass(az_deg: float) -> str:
    return COMPASS[int(((az_deg % 360) + 22.5) // 45) % 8]


def angular_sep(alt1, az1, alt2, az2):
    """Great-circle separation in radians between arrays of alt/az (radians)."""
    s = (np.sin(alt1) * np.sin(alt2)
         + np.cos(alt1) * np.cos(alt2) * np.cos(az1 - az2))
    return np.arccos(np.clip(s, -1.0, 1.0))


def geometry(ts, observer, obs_at, sat, body, body_radius_km, jd):
    t = ts.tt_jd(jd)
    s_alt, s_az, s_dist = (sat - observer).at(t).altaz()
    b_alt, b_az, b_dist = obs_at.at(t).observe(body).apparent().altaz()
    sep = angular_sep(s_alt.radians, s_az.radians, b_alt.radians, b_az.radians)
    radius = np.arcsin(body_radius_km / b_dist.km)
    return sep, radius, s_alt, s_az, s_dist, b_alt, b_az


def find_passes(ts, sat, observer, t0, t1, min_alt):
    """Yield (rise_jd, set_jd) for passes above min_alt, clipped to [t0, t1]."""
    times, events = sat.find_events(observer, t0, t1, altitude_degrees=min_alt)
    start = None
    if len(events) and events[0] != 0:
        start = t0.tt                      # pass already in progress at t0
    for t, e in zip(times, events):
        if e == 0:
            start = t.tt
        elif e == 2 and start is not None:
            yield start, t.tt
            start = None
    if start is not None:
        yield start, t1.tt


def scan(args, lat, lon, elev, tles, now_utc):
    ts = load.timescale()
    eph = load_ephemeris()
    earth, sun, moon = eph["earth"], eph["sun"], eph["moon"]
    observer = wgs84.latlon(lat, lon, elevation_m=elev)
    obs_at = earth + observer
    tz = ZoneInfo(args.tz)

    t0 = ts.from_datetime(now_utc)
    t1 = ts.from_datetime(now_utc + timedelta(days=args.days))
    bodies = [(n, b, r) for n, b, r in (("Sun", sun, SUN_RADIUS_KM),
                                         ("Moon", moon, MOON_RADIUS_KM))
              if n.lower() in args.bodies]
    found: list[Event] = []

    for catnr, label in STATIONS.items():
        name, l1, l2 = tles[catnr]
        sat = EarthSatellite(l1, l2, label, ts)
        tle_age = t0 - sat.epoch

        # A station inside the disc sits at the body's altitude, so a pass
        # only needs to clear (min_alt - 1) for the body to clear min_alt.
        for rise, set_ in find_passes(ts, sat, observer, t0, t1, args.min_alt - 1):
            mid = ts.tt_jd((rise + set_) / 2)
            for body_name, body, body_r in bodies:
                b_alt_mid = obs_at.at(mid).observe(body).apparent().altaz()[0].degrees
                if b_alt_mid < args.min_alt - 4:   # passes last < 10 min
                    continue

                n = max(3, int((set_ - rise) * 86400 / COARSE_STEP_S) + 1)
                jd = np.linspace(rise, set_, n)
                sep, *_ = geometry(ts, observer, obs_at, sat, body, body_r, jd)
                k = int(np.argmin(sep))

                lo = jd[max(k - 1, 0)]
                hi = jd[min(k + 1, n - 1)]
                m = max(3, int((hi - lo) * 86400 / FINE_STEP_S) + 1)
                jd_f = np.linspace(lo, hi, m)
                sep_f, rad_f, s_alt, s_az, s_dist, b_alt, b_az = geometry(
                    ts, observer, obs_at, sat, body, body_r, jd_f)
                j = int(np.argmin(sep_f))

                min_sep, radius = float(sep_f[j]), float(rad_f[j])
                if b_alt.degrees[j] < args.min_alt:
                    continue

                days_out = (jd_f[j] - t0.tt)
                if min_sep <= radius:
                    kind = "transit"
                elif (days_out >= args.possible_after_days
                      and min_sep <= radius * (1 + args.possible_margin)):
                    kind = "possible"
                else:
                    continue

                inside = sep_f <= rad_f
                duration = float(inside.sum()) * FINE_STEP_S
                if inside[0] or inside[-1]:
                    # Crossing runs past the refinement window; widen and recount.
                    jd_w = np.linspace(jd_f[j] - 3 / 86400, jd_f[j] + 3 / 86400,
                                       int(6 / FINE_STEP_S) + 1)
                    sep_w, rad_w, *_ = geometry(ts, observer, obs_at, sat, body, body_r, jd_w)
                    duration = float((sep_w <= rad_w).sum()) * FINE_STEP_S

                rng = float(s_dist.km[j])
                sin_alt = max(math.sin(s_alt.radians[j]), 0.05)
                when = ts.tt_jd(jd_f[j]).utc_datetime()
                found.append(Event(
                    station=label,
                    body=body_name,
                    kind=kind,
                    time_utc=when.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4] + "Z",
                    time_local=when.astimezone(tz).strftime("%a %b %d %Y %I:%M:%S %p %Z"),
                    days_out=round(days_out, 2),
                    min_sep_arcmin=round(math.degrees(min_sep) * 60, 2),
                    disc_radius_arcmin=round(math.degrees(radius) * 60, 2),
                    duration_s=round(duration, 2),
                    station_alt_deg=round(float(s_alt.degrees[j]), 1),
                    station_az_deg=round(float(s_az.degrees[j]), 1),
                    station_range_km=round(rng),
                    body_alt_deg=round(float(b_alt.degrees[j]), 1),
                    body_az_deg=round(float(b_az.degrees[j]), 1),
                    direction=compass(float(b_az.degrees[j])),
                    low=bool(b_alt.degrees[j] < args.low_alt),
                    offset_km=round(rng * min_sep / sin_alt, 2),
                    path_halfwidth_km=round(rng * radius / sin_alt, 2),
                    tle_age_days=round(float(tle_age), 2),
                ))
    found.sort(key=lambda e: e.time_utc)
    return found


# ---------------------------------------------------------------- output

def print_text(events, tle_source, args, now_utc, tle_age=None):
    print(f"Transit scan {now_utc:%Y-%m-%d %H:%MZ}, next {args.days} days, "
          f"body alt >= {args.min_alt} deg, TLEs: {tle_source}")
    if tle_age is not None and tle_age > 3:
        print(f"WARNING: TLEs are {tle_age:.1f} days old; predictions will drift.")
    if not events:
        print("No transits or possible transits found.")
        return
    for e in events:
        tag = "TRANSIT " if e.kind == "transit" else "POSSIBLE"
        line = (f"{tag} {e.station} across the {e.body}: {e.time_local} "
                f"({e.days_out} d out), {e.body} {e.body_alt_deg} deg up in the "
                f"{e.direction}{' (LOW, may be blocked)' if e.low else ''}, "
                f"station range {e.station_range_km} km")
        if e.kind == "transit":
            line += (f", {e.duration_s} s in disc, site {e.offset_km} km off "
                     f"centerline (path half-width {e.path_halfwidth_km} km)")
        else:
            line += (f", misses by {e.min_sep_arcmin - e.disc_radius_arcmin:.1f}' "
                     f"(~{e.offset_km - e.path_halfwidth_km:.1f} km off the path)")
        print(line)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lat", type=float, help="site latitude, degrees north")
    p.add_argument("--lon", type=float, help="site longitude, degrees east (west is negative)")
    p.add_argument("--elev", type=float, default=None, help="site elevation, meters")
    p.add_argument("--location-file", default=str(HERE / "location.json"))
    p.add_argument("--tle", help="TLE file to use instead of data/stations.tle")
    p.add_argument("--days", type=float, default=14)
    p.add_argument("--min-alt", type=float, default=10.0,
                   help="skip events with the Sun/Moon below this, degrees (default 10)")
    p.add_argument("--low-alt", type=float, default=20.0,
                   help="mark events below this as low/hard to see, degrees (default 20)")
    p.add_argument("--bodies", nargs="+", default=["sun", "moon"],
                   choices=["sun", "moon"])
    p.add_argument("--possible-after-days", type=float, default=3.0,
                   help="only flag near misses this many days out or more")
    p.add_argument("--possible-margin", type=float, default=1.0,
                   help="near-miss window, in disc radii beyond the limb")
    p.add_argument("--tz", default="America/Detroit")
    p.add_argument("--start", help="scan start, ISO UTC (default: now)")
    p.add_argument("--json", action="store_true", help="print JSON instead of text")
    args = p.parse_args(argv)

    lat, lon, elev = load_location(args)
    tles, tle_source = load_tles(args)
    now_utc = (datetime.fromisoformat(args.start.replace("Z", "+00:00"))
               if args.start else datetime.now(timezone.utc))
    if now_utc.tzinfo is None:
        now_utc = now_utc.replace(tzinfo=timezone.utc)

    events = scan(args, lat, lon, elev, tles, now_utc)
    if args.json:
        print(json.dumps({
            "scanned_at_utc": now_utc.isoformat(),
            "days": args.days,
            "tle_source": tle_source,
            "events": [asdict(e) for e in events],
        }, indent=2))
    else:
        ts = load.timescale()
        ages = [ts.from_datetime(now_utc) - EarthSatellite(l1, l2, "", ts).epoch
                for c, (_, l1, l2) in tles.items() if c in STATIONS]
        print_text(events, tle_source, args, now_utc, max(ages) if ages else None)


if __name__ == "__main__":
    main()
