# OpenGSIM-Bot

An automated workflow for the estimation and AI-powered interpretation of OpenQuake’s Ground-Motion Intensity Measures via a Telegram Bot.

This repository contains the source code and configuration files for the OpenGSIM-Bot, as presented in the manuscript  
**"Development of an Open-Access Telegram Bot for Automated Estimation and AI-Based Interpretation of OpenQuake Ground-Motion Intensity Measures"**  
submitted to *Computers & Geosciences*.

A live, interactive instance of the bot can be accessed on Telegram at:  
[https://t.me/MyPGAEstimator_bot](https://t.me/MyPGAEstimator_bot)

---

## Architecture

The system is a microservices-based application orchestrated by an **n8n** workflow.  
It is composed of four primary layers:

1. **Interface layer** – Telegram Bot API  
2. **Orchestration layer** – n8n  
3. **Computation layer** – Python wrapper for the OpenQuake engine (Docker container)  
4. **Intelligence layer** – Google’s Gemini API  

*(A technical diagram illustrating this architecture is available in the accompanying scientific paper.)*

---

## Deployment and Usage

The following four steps will guide you through the process of deploying and testing the application on your local machine.

### 1. Prerequisites

- Docker and Docker Compose installed on your system  
- A Telegram Bot token obtained from **@BotFather**  
- API keys for **Google Sheets** and **Google Gemini**

---

### 2. Configuration

Before launching, configure your environment variables.

#### Create a `.env` file

Copy the example file `.env.example` to a new file named `.env`:

```bash
cp .env.example .env
```

#### Edit the `.env` file

Open the newly created `.env` file with a text editor and fill in the required credentials for Telegram, Google Sheets, and Gemini.

#### Set up Google Sheet

1. Use the `gmm_catalog.csv_for_main_workflow.xlsx` file provided in the `/data` directory as a template for your own Google Sheet.  
2. Share your new sheet with the `client_email` found in your Google Service Account credentials.  
3. Update the “Sheet ID” in the relevant nodes of the **n8n** workflow (`Telegram_GMM_ AI agents_google_sheet.json`) to point to your copy.

---

### 3. Launching the Application

This is the main step to start all necessary services (**n8n** and the **OpenQuake** wrapper).

From the root directory of the project, run:

```bash
docker-compose up -d
```

This command will build the necessary Docker images and start both the n8n and OpenQuake services in the background (`-d` for detached mode).  
Please wait a minute or two for the services to initialize completely.  
You can check the status of the containers by running:

```bash
docker ps
```

---

### 4. Quick Test

Once the services are running, you can verify that the scientific computation layer is working correctly.  
**Important:** You must perform **“Launching the Application” (Step 3)** before running this test.

The provided test script will execute a single GMM calculation for a predefined scenario.  
From the root directory of the project, run:

```bash
# Make the script executable (only needs to be done once)
chmod +x quick_test.sh

# Run the test
./quick_test.sh
```

A successful test will print a large JSON object to the console with the key `"success": true`.  
This confirms that the scientific core of the application is functioning correctly.

---

## Citation

If you use this work in your research, please cite our accompanying paper:

> **Derras, B. (Year).** Development of an Open-Access Telegram Bot for Automated Estimation and AI-Based Interpretation of OpenQuake Ground-Motion Intensity Measures.  
> *Computers & Geosciences.* (Details to be updated upon publication).

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.
