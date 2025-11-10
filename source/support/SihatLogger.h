#pragma once

#include<string>
#include<vector>
#include<mutex>

class SihatLogger 
{
public:
	void log(const std::string& message)
	{
		std::lock_guard<std::mutex> lock(m_mutex);
		permanentMessages.push_back(message);
	}

	void logTemp(const std::string& message)
	{
		std::lock_guard<std::mutex> lock(m_mutex);
		temporaryMessage = message;
	}

	std::mutex& getMutex() { return m_mutex; }
	const std::string& getTemporaryMessage() const { return temporaryMessage; }
	const std::vector<std::string>& getPermanentMessages() const { return permanentMessages; }




private:
	std::mutex m_mutex;
	std::string temporaryMessage;
	std::vector<std::string> permanentMessages;
};