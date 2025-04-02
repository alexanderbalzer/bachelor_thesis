from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import time

# Set up Selenium WebDriver (make sure you have the appropriate driver for your browser)
driver = webdriver.Chrome(executable_path='/path/to/chromedriver')  # Adjust to your path

# Open the MitoFates page
driver.get("https://mitf.cbrc.pj.aist.go.jp/MitoFates/cgi-bin/top.cgi")

# Find the protein sequence input field (you can inspect the HTML for the correct ID or name)
sequence_input = driver.find_element(By.NAME, "sequence")  # Adjust according to the form's name attribute

# Enter the sequence
sequence_input.send_keys("YOUR_PROTEIN_SEQUENCE")  # Replace with your sequence

# Submit the form
sequence_input.send_keys(Keys.RETURN)  # This simulates pressing the "Enter" key

# Wait for the results to load (adjust the time as needed)
time.sleep(5)

# Extract the results (inspect the HTML for where the results are displayed)
results = driver.find_elements(By.CLASS_NAME, "result_class")  # Example, adjust as needed
for result in results:
    print(result.text)

# Close the browser
driver.quit()
