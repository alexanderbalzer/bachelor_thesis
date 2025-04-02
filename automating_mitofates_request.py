from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time


names = ["human"]
name = names[0]

# Initialize the Firefox driver
driver = webdriver.Firefox()

try:
    # The file you want to upload
    file_to_upload = f"/home/abalzer/Documents/github_clone/bachelor_thesis/output files/filtered_proteins_by_GO_for_{name}.fasta"

    # Open the MitoFates website
    driver.get("https://mitf.cbrc.pj.aist.go.jp/MitoFates/cgi-bin/top.cgi")

    # Wait for the file input element to be present
    file_input = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, "queryFile"))
    )
    file_input.send_keys(file_to_upload)  # Upload the file

    # Select "metazoa" from the dropdown
    select_element = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, "organism"))
    )
    select = Select(select_element)
    select.select_by_visible_text("metazoa")

    # Find the submit button and click it
    submit_button = WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.NAME, "submit"))
    )
    submit_button.click()

    # Wait for the results to load
    results_element = WebDriverWait(driver, 30).until(
        EC.presence_of_element_located((By.ID, "results"))  # Replace "results" with the actual ID or locator
    )
    print(results_element.text)  # Print the results

finally:
    # Close the browser
    driver.quit()