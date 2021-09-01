#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 01:14:31 2021

@author: nicolaskylilis
"""
def send_email(script_name=""):
    
    import smtplib, ssl
    
    port = 465  # For SSL
    password = "nicolaspython27@server"
    
    # Create a secure SSL context
    context = ssl.create_default_context()
    
    with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
        server.login("nicolaspython27@gmail.com", password)
        
        # TODO: Send email here
        server.sendmail("nicolaspython27@gmail.com","kylilis.nicolas@ucy.ac.cy",script_name + ": JOB DONE!")
        
