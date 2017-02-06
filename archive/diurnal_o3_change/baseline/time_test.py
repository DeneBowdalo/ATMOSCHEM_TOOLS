import pytz, datetime
local = pytz.timezone ("Atlantic/Cape_Verde")
naive = datetime.datetime.strptime ("2001-2-3 10:11:12", "%Y-%m-%d %H:%M:%S")
local_dt = local.localize(naive, is_dst=None)
utc_dt = local_dt.astimezone (pytz.utc)

a = utc_dt.strftime ("%Y-%m-%d %H:%M:%S")

print a
