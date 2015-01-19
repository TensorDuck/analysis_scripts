"""
Function for appending the log from recipe_manager, so that it may be called separately.
"""
import os
import time

class log_function(object):
    def __init__(self, cwd):
        self.path = cwd

    def append_log(self,sub,string,subdir=False):
        if subdir:
            open('%s/%s/%s.log' % (self.path,sub,sub),'a').write("%s %s\n" % (self.append_time_label(),string))
        else:
            open('%s/%s/modelbuilder.log' % (self.path,sub),'a').write("%s %s\n" % (self.append_time_label(),string))

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string
